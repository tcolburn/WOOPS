#!/usr/bin/env python3
"""Combine scalar Kirkwood leg time series and apply optional equilibration trimming.

Designed for the dedicated scalar Kirkwood WOOPS executable.  The Fortran
program writes one compact CSV-formatted .dat file per region and per DCD leg.
This script:

  * discovers/sorts three-digit leg directories,
  * reads trajectory_leg_lengths.dat for per-leg n_frames and trajectory_ns,
  * reconstructs cumulative simulation time across legs,
  * writes one complete Parquet time series per Kirkwood region,
  * keeps zero-occupancy frames in the time series,
  * excludes valid=0 frames from statistical averages,
  * applies --discard-ns only to the reported summary statistics.

The raw Fortran outputs are never modified.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import polars as pl


REGION_FILES = {
    "interface_7A": "kirkwood_interface.dat",
    "PRO16_7A": "kirkwood_PRO16_7A.dat",
    "PRO16_10A": "kirkwood_PRO16_10A.dat",
}

EXPECTED_COLUMNS = [
    "frame",
    "valid",
    "N",
    "normal_coherence",
    "g_k",
    "g_perp",
    "g_parallel",
    "decomposition_residual",
]


@dataclass(frozen=True)
class LegMeta:
    leg: str
    n_frames: int
    trajectory_ns: float
    frame_dt_ns: float
    frame_offset: int
    time_offset_ns: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine dedicated scalar Kirkwood WOOPS leg outputs and calculate "
            "trajectory averages after an optional initial-time discard."
        )
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Simulation-specific Kirkwood directory. Default: current directory.",
    )
    parser.add_argument(
        "--length-table",
        type=Path,
        default=None,
        help=(
            "trajectory_leg_lengths.dat. Default: search ROOT and its ancestors "
            "for the nearest file with this name."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: ROOT/kirkwood_processed.",
    )
    parser.add_argument(
        "--discard-ns",
        type=float,
        default=0.0,
        help=(
            "Discard this initial amount of cumulative simulation time from "
            "summary averages. Full combined Parquet files are never trimmed. "
            "Default: 0."
        ),
    )
    parser.add_argument(
        "--legs",
        nargs="+",
        default=None,
        help="Optional explicit three-digit leg list, e.g. --legs 001 002 003.",
    )
    parser.add_argument(
        "--skip-incomplete",
        action="store_true",
        help=(
            "Skip numbered leg directories missing one or more Kirkwood region "
            "files. Default behavior is to stop so missing production data are "
            "not silently ignored."
        ),
    )
    parser.add_argument(
        "--allow-thinned",
        action="store_true",
        help=(
            "Allow time-series files with fewer rows than trajectory_leg_lengths.dat. "
            "Use only when the Fortran timeSeriesStride was intentionally > 1."
        ),
    )
    return parser.parse_args()


def normalize_simulation(value: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError("Empty simulation identifier.")
    if text.lower() in {"wildtype", "wt"}:
        return "Wildtype"
    if re.fullmatch(r"\d+", text):
        return f"{int(text):03d}"
    return text


def normalize_leg(value: str) -> str:
    text = str(value).strip()
    if not re.fullmatch(r"\d+", text):
        raise ValueError(f"Invalid numeric leg identifier: {value!r}")
    number = int(text)
    if number < 0 or number > 999:
        raise ValueError(f"Leg identifier out of range: {value!r}")
    return f"{number:03d}"


def find_length_table(root: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        path = explicit.expanduser()
        if not path.is_absolute():
            path = root / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Trajectory-length table not found: {path}")
        return path

    for directory in [root, *root.parents]:
        candidate = directory / "trajectory_leg_lengths.dat"
        if candidate.is_file():
            return candidate.resolve()

    raise FileNotFoundError(
        "Could not find trajectory_leg_lengths.dat in ROOT or any ancestor. "
        "Pass it explicitly with --length-table."
    )


def load_leg_metadata(path: Path, simulation_id: str) -> dict[str, LegMeta]:
    required = {"directory", "subdirectory", "trajectory_ns", "n_frames"}
    rows: list[tuple[str, int, float]] = []

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")
        missing = required.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                f"{path} is missing required columns: {sorted(missing)}; "
                f"found {reader.fieldnames}"
            )

        for line_number, row in enumerate(reader, start=2):
            raw_directory = (row.get("directory") or "").strip()
            if not raw_directory:
                continue
            try:
                directory = normalize_simulation(raw_directory)
            except ValueError:
                continue
            if directory != simulation_id:
                continue

            raw_leg = (row.get("subdirectory") or "").strip()
            if not raw_leg:
                raise ValueError(
                    f"Missing subdirectory for simulation {simulation_id} "
                    f"at {path}:{line_number}."
                )
            leg = normalize_leg(raw_leg)

            try:
                n_frames_float = float((row.get("n_frames") or "").strip())
                trajectory_ns = float((row.get("trajectory_ns") or "").strip())
            except ValueError as exc:
                raise ValueError(
                    f"Invalid n_frames/trajectory_ns at {path}:{line_number}."
                ) from exc

            n_frames = int(round(n_frames_float))
            if not math.isclose(n_frames_float, n_frames, rel_tol=0.0, abs_tol=0.0):
                raise ValueError(
                    f"n_frames must be integer-valued at {path}:{line_number}: "
                    f"{n_frames_float}"
                )
            if n_frames <= 0 or trajectory_ns <= 0.0:
                raise ValueError(
                    f"n_frames and trajectory_ns must be positive at "
                    f"{path}:{line_number}."
                )

            rows.append((leg, n_frames, trajectory_ns))

    if not rows:
        raise ValueError(
            f"No rows for simulation {simulation_id} were found in {path}."
        )

    duplicate_legs = sorted(
        {leg for leg, _, _ in rows if sum(r[0] == leg for r in rows) > 1},
        key=int,
    )
    if duplicate_legs:
        raise ValueError(
            f"Duplicate trajectory-length rows for simulation {simulation_id}: "
            f"{duplicate_legs}"
        )

    metadata: dict[str, LegMeta] = {}
    frame_offset = 0
    time_offset_ns = 0.0
    for leg, n_frames, trajectory_ns in sorted(rows, key=lambda item: int(item[0])):
        metadata[leg] = LegMeta(
            leg=leg,
            n_frames=n_frames,
            trajectory_ns=trajectory_ns,
            frame_dt_ns=trajectory_ns / n_frames,
            frame_offset=frame_offset,
            time_offset_ns=time_offset_ns,
        )
        frame_offset += n_frames
        time_offset_ns += trajectory_ns

    return metadata


def discover_legs(root: Path, requested: list[str] | None) -> list[str]:
    available = sorted(
        [p.name for p in root.iterdir() if p.is_dir() and re.fullmatch(r"\d{3}", p.name)],
        key=int,
    )
    if not available:
        raise RuntimeError(f"No three-digit leg directories found in {root}.")

    if requested is None:
        return available

    legs = sorted({normalize_leg(value) for value in requested}, key=int)
    missing = [leg for leg in legs if leg not in available]
    if missing:
        raise FileNotFoundError(f"Requested leg directories do not exist: {missing}")
    return legs


def read_region_file(path: Path) -> pl.DataFrame:
    dataframe = pl.read_csv(
        path,
        comment_prefix="#",
        schema_overrides={
            "frame": pl.Int64,
            "valid": pl.Int8,
            "N": pl.Int64,
            "normal_coherence": pl.Float64,
            "g_k": pl.Float64,
            "g_perp": pl.Float64,
            "g_parallel": pl.Float64,
            "decomposition_residual": pl.Float64,
        },
    )

    if dataframe.columns != EXPECTED_COLUMNS:
        raise ValueError(
            f"Unexpected columns in {path}: {dataframe.columns}; "
            f"expected {EXPECTED_COLUMNS}."
        )
    if dataframe.height == 0:
        raise ValueError(f"No data rows found in {path}.")

    bad_valid = dataframe.filter(~pl.col("valid").is_in([0, 1]))
    if bad_valid.height:
        raise ValueError(f"{path} contains valid values other than 0 or 1.")

    bad_occupancy = dataframe.filter(
        ((pl.col("valid") == 0) & (pl.col("N") != 0))
        | ((pl.col("valid") == 1) & (pl.col("N") <= 0))
    )
    if bad_occupancy.height:
        raise ValueError(
            f"{path} has inconsistent valid/N values; first bad row: "
            f"{bad_occupancy.head(1).to_dicts()[0]}"
        )

    return dataframe


def validate_frames(
    dataframe: pl.DataFrame,
    path: Path,
    meta: LegMeta,
    allow_thinned: bool,
) -> None:
    frames = dataframe.get_column("frame")
    first_frame = int(frames[0])
    last_frame = int(frames[-1])

    if first_frame < 1 or last_frame > meta.n_frames:
        raise ValueError(
            f"{path}: frame range {first_frame}..{last_frame} lies outside "
            f"trajectory_leg_lengths.dat range 1..{meta.n_frames}."
        )

    if frames.n_unique() != dataframe.height:
        raise ValueError(f"{path}: duplicate frame numbers found.")

    if not bool(frames.is_sorted()):
        raise ValueError(f"{path}: frame numbers are not monotonically increasing.")

    if not allow_thinned:
        if dataframe.height != meta.n_frames:
            raise ValueError(
                f"{path}: contains {dataframe.height} rows but the trajectory table "
                f"reports {meta.n_frames} frames. This looks truncated or thinned. "
                "Use --allow-thinned only if timeSeriesStride > 1 was intentional."
            )
        expected = pl.Series("frame", range(1, meta.n_frames + 1), dtype=pl.Int64)
        if not frames.equals(expected):
            raise ValueError(
                f"{path}: full time series does not contain exactly frames "
                f"1..{meta.n_frames}."
            )


def cutoff_label(discard_ns: float) -> str:
    if math.isclose(discard_ns, round(discard_ns), rel_tol=0.0, abs_tol=1.0e-12):
        return f"{int(round(discard_ns))}ns"
    return f"{discard_ns:g}ns".replace(".", "p")


def summarize_region(
    dataframe: pl.DataFrame,
    region: str,
    simulation_id: str,
    discard_ns: float,
) -> dict[str, object]:
    after_cut = dataframe.filter(pl.col("time_ns") > discard_ns)
    valid = after_cut.filter(pl.col("valid") == 1)

    if after_cut.height == 0:
        raise ValueError(
            f"Discard cutoff {discard_ns:g} ns removes every frame in region {region}."
        )
    if valid.height == 0:
        raise ValueError(
            f"No occupied/valid frames remain after {discard_ns:g} ns in region {region}."
        )

    values = valid.select(
        pl.col("N").mean().alias("mean_N"),
        pl.col("normal_coherence").mean().alias("mean_normal_coherence"),
        pl.col("g_k").mean().alias("g_k"),
        pl.col("g_perp").mean().alias("g_perp"),
        pl.col("g_parallel").mean().alias("g_parallel"),
        pl.col("decomposition_residual").mean().alias("decomposition_residual"),
    ).row(0, named=True)

    g_k = float(values["g_k"])
    ratio_perp = float(values["g_perp"]) / g_k if abs(g_k) > 1.0e-15 else math.nan
    ratio_parallel = (
        float(values["g_parallel"]) / g_k if abs(g_k) > 1.0e-15 else math.nan
    )

    first_time = float(after_cut.get_column("time_ns").min())
    last_time = float(after_cut.get_column("time_ns").max())

    return {
        "simulation": simulation_id,
        "region": region,
        "discard_ns": discard_ns,
        "frames_after_discard": after_cut.height,
        "valid_frames": valid.height,
        "zero_occupancy_frames": after_cut.height - valid.height,
        "valid_fraction": valid.height / after_cut.height,
        "first_included_time_ns": first_time,
        "last_included_time_ns": last_time,
        "mean_N": float(values["mean_N"]),
        "mean_normal_coherence": float(values["mean_normal_coherence"]),
        "g_k": g_k,
        "g_perp": float(values["g_perp"]),
        "g_parallel": float(values["g_parallel"]),
        "decomposition_residual": float(values["decomposition_residual"]),
        "g_perp_over_g_k": ratio_perp,
        "g_parallel_over_g_k": ratio_parallel,
        "g_perp_ratio_minus_one_third": ratio_perp - 1.0 / 3.0,
        "g_parallel_ratio_minus_two_thirds": ratio_parallel - 2.0 / 3.0,
    }


def main() -> None:
    args = parse_args()
    if args.discard_ns < 0.0:
        raise ValueError("--discard-ns must be >= 0.")

    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Kirkwood root directory not found: {root}")

    simulation_id = normalize_simulation(root.name)
    length_table = find_length_table(root, args.length_table)
    metadata = load_leg_metadata(length_table, simulation_id)

    output_dir = (
        args.output_dir.expanduser()
        if args.output_dir is not None
        else root / "kirkwood_processed"
    )
    if not output_dir.is_absolute():
        output_dir = root / output_dir
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    candidate_legs = discover_legs(root, args.legs)
    selected_legs: list[str] = []
    skipped: list[tuple[str, list[str]]] = []

    for leg in candidate_legs:
        if leg not in metadata:
            raise ValueError(
                f"Leg {leg} exists under {root} but has no trajectory-length row "
                f"for simulation {simulation_id} in {length_table}."
            )
        missing = [
            filename
            for filename in REGION_FILES.values()
            if not (root / leg / filename).is_file()
        ]
        if missing:
            if args.skip_incomplete:
                skipped.append((leg, missing))
                continue
            raise FileNotFoundError(
                f"Leg {leg} is incomplete; missing {missing}. "
                "Use --skip-incomplete only if omission is intentional."
            )
        selected_legs.append(leg)

    if not selected_legs:
        raise RuntimeError("No complete Kirkwood leg directories selected.")

    if skipped:
        print("Skipping incomplete legs:")
        for leg, missing in skipped:
            print(f"  {leg}: {', '.join(missing)}")

    print(f"Simulation: {simulation_id}")
    print(f"Root: {root}")
    print(f"Trajectory table: {length_table}")
    print(f"Legs: {' '.join(selected_legs)}")
    print(f"Summary discard: {args.discard_ns:g} ns")

    region_frames: dict[str, list[pl.DataFrame]] = {region: [] for region in REGION_FILES}
    leg_manifest_rows: list[dict[str, object]] = []

    for leg in selected_legs:
        meta = metadata[leg]
        reference_frames: pl.Series | None = None

        for region, filename in REGION_FILES.items():
            path = root / leg / filename
            dataframe = read_region_file(path)
            validate_frames(dataframe, path, meta, args.allow_thinned)

            frames = dataframe.get_column("frame")
            if reference_frames is None:
                reference_frames = frames
            elif not frames.equals(reference_frames):
                raise ValueError(
                    f"Frame indices differ among region files in leg {leg}; "
                    f"problem detected in {path}."
                )

            dataframe = dataframe.rename({"frame": "local_frame"}).with_columns(
                pl.lit(leg).alias("leg"),
                (pl.col("local_frame") + meta.frame_offset).alias("global_frame"),
                (pl.col("local_frame") * meta.frame_dt_ns).alias("local_time_ns"),
                (
                    pl.lit(meta.time_offset_ns)
                    + pl.col("local_frame") * meta.frame_dt_ns
                ).alias("time_ns"),
            ).select(
                "leg",
                "local_frame",
                "global_frame",
                "local_time_ns",
                "time_ns",
                "valid",
                "N",
                "normal_coherence",
                "g_k",
                "g_perp",
                "g_parallel",
                "decomposition_residual",
            )

            region_frames[region].append(dataframe)

        leg_manifest_rows.append(
            {
                "simulation": simulation_id,
                "leg": leg,
                "n_frames_table": meta.n_frames,
                "trajectory_ns": meta.trajectory_ns,
                "frame_dt_ns": meta.frame_dt_ns,
                "frame_offset": meta.frame_offset,
                "time_offset_ns": meta.time_offset_ns,
            }
        )

    summary_rows: list[dict[str, object]] = []
    for region, pieces in region_frames.items():
        combined = pl.concat(pieces, how="vertical").sort("global_frame")
        parquet_path = output_dir / f"kirkwood_{region}_timeseries.parquet"
        combined.write_parquet(parquet_path, compression="zstd")
        summary_rows.append(
            summarize_region(combined, region, simulation_id, args.discard_ns)
        )
        print(f"Wrote {parquet_path} ({combined.height} rows)")

    manifest = pl.DataFrame(leg_manifest_rows)
    manifest_path = output_dir / "kirkwood_leg_metadata.tsv"
    manifest.write_csv(manifest_path, separator="\t")

    summary = pl.DataFrame(summary_rows)
    summary_path = output_dir / f"kirkwood_summary_after_{cutoff_label(args.discard_ns)}.tsv"
    summary.write_csv(summary_path, separator="\t")

    print(f"Wrote {manifest_path}")
    print(f"Wrote {summary_path}")
    print("Raw per-leg Kirkwood files were not modified.")


if __name__ == "__main__":
    main()
