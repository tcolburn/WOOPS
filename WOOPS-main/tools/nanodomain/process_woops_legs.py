#!/usr/bin/env python3
"""
Collect and average WOOPS p1p2.dat output from independent trajectory legs.

Run from a simulation-specific WOOPS directory such as:

    /path/to/WOOPS/012
    /path/to/WOOPS/Wildtype

The simulation directory may be a positive three-digit mutation directory or
the special Wildtype/WT directory. The script auto-detects positive three-digit
leg directories, validates each
p1p2.dat against trajectory_leg_lengths.dat, calculates frame-weighted p1,
p2, and p21 averages, writes mapped PSFs, and writes a VMD render config.

Nothing inside the leg directories is renamed or modified.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

P1P2_COLUMNS = [
    "woops_frames",
    "atom_id",
    "resid",
    "atom_name",
    "resname",
    "original_charge",
    "p1",
    "p2",
    "p21",
]

METADATA_COLUMNS = [
    "atom_id",
    "resid",
    "atom_name",
    "resname",
    "original_charge",
]

VALUE_COLUMNS = ["p1", "p2", "p21"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Frame-weight and collect WOOPS p1p2.dat files."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Simulation-specific WOOPS directory. Default: current directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Default: ROOT/woops_processed.",
    )
    parser.add_argument(
        "--length-table",
        type=Path,
        default=None,
        help=(
            "Tab-separated trajectory_leg_lengths.dat. "
            "Default: ROOT/../../trajectory_leg_lengths.dat."
        ),
    )
    parser.add_argument(
        "--psf",
        type=Path,
        default=None,
        help="Original root-level PSF. Default: auto-detect the sole ROOT/*.psf.",
    )
    parser.add_argument(
        "--pdb",
        type=Path,
        default=None,
        help=(
            "Coordinate PDB. Default: ROOT/1YLB_Small_waters.pdb when present; "
            "otherwise auto-detect a root-level non-beta PDB."
        ),
    )
    parser.add_argument(
        "--legs",
        nargs="+",
        default=None,
        help="Optional explicit legs, e.g. --legs 001 002 003.",
    )
    parser.add_argument(
        "--p1p2-name",
        default="p1p2.dat",
        help="Per-leg WOOPS filename. Default: p1p2.dat.",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help=(
            "Skip numbered leg directories lacking p1p2.dat. "
            "By default, their absence is treated as an error."
        ),
    )
    return parser.parse_args()


def normalize_three_digit(value: object, label: str) -> str:
    text = str(value).strip()
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    if not text.isdigit():
        raise ValueError(f"{label} must be numeric; received {value!r}.")
    number = int(text)
    if number <= 0:
        raise ValueError(f"{label} must be positive; received {value!r}.")
    return f"{number:03d}"


def normalize_simulation_id(value: object, label: str) -> str:
    """Normalize mutation directories and the special wildtype dataset.

    Numeric simulation IDs retain the existing three-digit convention. Both
    'Wildtype' and 'WT' are normalized to 'WT', matching the directory label in
    trajectory_leg_lengths.dat and the generated output metadata.
    """
    text = str(value).strip()
    if text.lower() in {"wildtype", "wt"}:
        return "WT"
    return normalize_three_digit(value, label)


def parse_boolean(value: object, label: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "t", "1", "yes", "y"}:
        return True
    if text in {"false", "f", "0", "no", "n"}:
        return False
    raise ValueError(f"{label} must be boolean; received {value!r}.")


def numbered_subdirectories(root: Path) -> list[Path]:
    return sorted(
        (
            path
            for path in root.iterdir()
            if path.is_dir()
            and len(path.name) == 3
            and path.name.isdigit()
            and int(path.name) > 0
        ),
        key=lambda path: int(path.name),
    )


def resolve_single_file(
    explicit: Path | None,
    root: Path,
    pattern: str,
    description: str,
    exclude_names: Iterable[str] = (),
    preferred_name: str | None = None,
) -> Path:
    if explicit is not None:
        path = explicit.expanduser()
        if not path.is_absolute():
            path = root / path
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(f"{description} not found: {path}")
        return path

    if preferred_name is not None:
        preferred = root / preferred_name
        if preferred.is_file():
            return preferred.resolve()

    excluded = {name.lower() for name in exclude_names}
    candidates = sorted(
        path.resolve()
        for path in root.glob(pattern)
        if path.is_file() and path.name.lower() not in excluded
    )

    if len(candidates) != 1:
        names = ", ".join(path.name for path in candidates) or "none"
        raise RuntimeError(
            f"Could not uniquely auto-detect {description}. "
            f"Candidates matching {pattern!r}: {names}. Pass the file explicitly."
        )
    return candidates[0]


def parse_p1p2(path: Path, leg: str) -> pd.DataFrame:
    records: list[dict[str, object]] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.lstrip().startswith("DATA"):
                continue

            tokens = line.split()
            try:
                colon_index = tokens.index(":")
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: expected a ':' token in WOOPS output."
                ) from exc

            payload = tokens[colon_index + 1 :]
            if len(payload) != 9:
                raise ValueError(
                    f"{path}:{line_number}: expected 9 fields after ':', "
                    f"found {len(payload)}: {payload}"
                )

            try:
                record = {
                    "woops_frames": int(payload[0]),
                    "atom_id": int(payload[1]),
                    "resid": payload[2].strip(),
                    "atom_name": payload[3].strip(),
                    "resname": payload[4].strip(),
                    "original_charge": float(payload[5]),
                    "p1": float(payload[6]),
                    "p2": float(payload[7]),
                    "p21": float(payload[8]),
                }
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{line_number}: could not parse numeric WOOPS fields."
                ) from exc

            records.append(record)

    if not records:
        raise ValueError(f"No 'DATA DI :' records found in {path}.")

    frame_counts = {int(record["woops_frames"]) for record in records}
    if len(frame_counts) != 1:
        raise ValueError(
            f"{path} contains inconsistent WOOPS frame counts: {sorted(frame_counts)}"
        )

    dataframe = pd.DataFrame.from_records(records, columns=P1P2_COLUMNS)
    dataframe.insert(0, "leg", leg)
    dataframe.insert(1, "source_file", str(path.resolve()))

    duplicate_ids = dataframe.loc[dataframe["atom_id"].duplicated(), "atom_id"].tolist()
    if duplicate_ids:
        raise ValueError(
            f"{path} contains duplicate atom IDs, including {duplicate_ids[:10]}."
        )

    if not np.isfinite(dataframe[VALUE_COLUMNS].to_numpy(dtype=float)).all():
        raise ValueError(f"{path} contains NaN or infinite p1/p2/p21 values.")

    tolerance = 1.0e-8
    bounds = {
        "p1": (-1.0 - tolerance, 1.0 + tolerance),
        "p2": (-0.5 - tolerance, 1.0 + tolerance),
        "p21": (-0.5 - tolerance, 0.5 + tolerance),
    }
    for column, (lower, upper) in bounds.items():
        bad = dataframe.loc[
            (dataframe[column] < lower) | (dataframe[column] > upper),
            ["atom_id", column],
        ]
        if not bad.empty:
            example = bad.head(5).to_dict(orient="records")
            raise ValueError(
                f"{path} contains {column} values outside [{lower}, {upper}]. "
                f"Examples: {example}"
            )

    return dataframe.sort_values("atom_id").reset_index(drop=True)


def load_length_table(path: Path, simulation_id: str) -> pd.DataFrame:
    """Load and validate only the rows needed for the current simulation.

    Filtering by simulation ID occurs before numeric conversion. This prevents
    incomplete or diagnostic rows belonging to unrelated simulations from
    causing a global ``NaN -> int`` conversion failure.
    """
    if not path.is_file():
        raise FileNotFoundError(
            f"Trajectory-length table not found: {path}\n"
            "Pass it explicitly with --length-table."
        )

    dataframe = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
    )

    required = {
        "directory",
        "subdirectory",
        "delta_mu_D",
        "trajectory_ns",
        "n_frames",
        "all_subdirs_have_dcd",
    }
    missing = required.difference(dataframe.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {sorted(missing)}. "
            f"Found: {list(dataframe.columns)}"
        )

    dataframe = dataframe.copy()
    dataframe["_table_line"] = np.arange(2, len(dataframe) + 2, dtype=np.int64)

    def normalize_directory_cell(value: object) -> str | None:
        if pd.isna(value) or not str(value).strip():
            return None
        return normalize_simulation_id(value, "directory")

    dataframe["directory"] = dataframe["directory"].map(normalize_directory_cell)

    # Parse only the current simulation. Other rows may legitimately be
    # incomplete while separate calculations are still in progress.
    dataframe = dataframe.loc[dataframe["directory"] == simulation_id].copy()
    if dataframe.empty:
        raise ValueError(
            f"No directory={simulation_id} rows found in {path}."
        )

    def bad_rows_message(column: str, mask: pd.Series) -> str:
        columns = [
            "_table_line",
            "directory",
            "subdirectory",
            "delta_mu_D",
            "trajectory_ns",
            "n_frames",
            "all_subdirs_have_dcd",
        ]
        rows = dataframe.loc[mask, columns].to_dict(orient="records")
        return f"Invalid {column} values in {path}: {rows}"

    missing_subdirectory = (
        dataframe["subdirectory"].isna()
        | dataframe["subdirectory"].astype(str).str.strip().eq("")
    )
    if missing_subdirectory.any():
        raise ValueError(bad_rows_message("subdirectory", missing_subdirectory))

    dataframe["subdirectory"] = dataframe["subdirectory"].map(
        lambda value: normalize_three_digit(value, "subdirectory")
    )

    numeric_columns = ["n_frames", "trajectory_ns", "delta_mu_D"]
    for column in numeric_columns:
        parsed = pd.to_numeric(dataframe[column], errors="coerce")
        bad = ~np.isfinite(parsed.to_numpy(dtype=float))
        if bad.any():
            raise ValueError(
                bad_rows_message(column, pd.Series(bad, index=dataframe.index))
            )
        dataframe[column] = parsed

    non_integer_frames = ~np.isclose(
        dataframe["n_frames"].to_numpy(dtype=float),
        np.rint(dataframe["n_frames"].to_numpy(dtype=float)),
        rtol=0.0,
        atol=0.0,
    )
    if non_integer_frames.any():
        raise ValueError(
            bad_rows_message(
                "n_frames (must be integer-valued)",
                pd.Series(non_integer_frames, index=dataframe.index),
            )
        )

    dataframe["n_frames"] = np.rint(dataframe["n_frames"]).astype(np.int64)

    try:
        dataframe["all_subdirs_have_dcd"] = dataframe[
            "all_subdirs_have_dcd"
        ].map(lambda value: parse_boolean(value, "all_subdirs_have_dcd"))
    except ValueError as exc:
        raise ValueError(
            f"{path}: invalid all_subdirs_have_dcd value among "
            f"directory={simulation_id} rows."
        ) from exc

    duplicates = dataframe.duplicated(["directory", "subdirectory"], keep=False)
    if duplicates.any():
        rows = dataframe.loc[
            duplicates, ["_table_line", "directory", "subdirectory"]
        ]
        raise ValueError(
            f"{path} contains duplicate directory/subdirectory rows: "
            f"{rows.to_dict(orient='records')}"
        )

    return dataframe.drop(columns="_table_line")


def validate_leg_metadata(reference: pd.DataFrame, current: pd.DataFrame, path: Path) -> None:
    if len(reference) != len(current):
        raise ValueError(
            f"{path} contains {len(current)} atoms; expected {len(reference)}."
        )

    left = reference[METADATA_COLUMNS].sort_values("atom_id").reset_index(drop=True)
    right = current[METADATA_COLUMNS].sort_values("atom_id").reset_index(drop=True)

    exact_columns = ["atom_id", "resid", "atom_name", "resname"]
    for column in exact_columns:
        mismatch = left[column].astype(str) != right[column].astype(str)
        if mismatch.any():
            index = int(np.flatnonzero(mismatch.to_numpy())[0])
            raise ValueError(
                f"{path} metadata mismatch in {column} at row {index}: "
                f"expected {left.loc[index, column]!r}, "
                f"found {right.loc[index, column]!r}."
            )

    charge_mismatch = ~np.isclose(
        left["original_charge"].to_numpy(dtype=float),
        right["original_charge"].to_numpy(dtype=float),
        rtol=0.0,
        atol=5.0e-8,
    )
    if charge_mismatch.any():
        index = int(np.flatnonzero(charge_mismatch)[0])
        raise ValueError(
            f"{path} original-charge mismatch for atom "
            f"{int(left.loc[index, 'atom_id'])}: "
            f"expected {left.loc[index, 'original_charge']}, "
            f"found {right.loc[index, 'original_charge']}."
        )


def read_psf(path: Path) -> tuple[list[str], dict[int, dict[str, object]], list[int]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    natom_header_index = None
    natom = None
    for index, line in enumerate(lines):
        if "!NATOM" in line:
            match = re.match(r"\s*(\d+)", line)
            if match is None:
                raise ValueError(f"Could not parse !NATOM count in {path}: {line!r}")
            natom_header_index = index
            natom = int(match.group(1))
            break

    if natom_header_index is None or natom is None:
        raise ValueError(f"No !NATOM section found in {path}.")

    records: dict[int, dict[str, object]] = {}
    atom_line_indices: list[int] = []
    cursor = natom_header_index + 1

    while cursor < len(lines) and len(atom_line_indices) < natom:
        line = lines[cursor]
        if line.strip():
            tokens = line.split()
            if len(tokens) < 8:
                raise ValueError(
                    f"Malformed PSF atom record in {path} at line {cursor + 1}: {line!r}"
                )
            try:
                atom_id = int(tokens[0])
                charge = float(tokens[6])
                mass = float(tokens[7])
            except ValueError as exc:
                raise ValueError(
                    f"Malformed numeric PSF atom fields in {path} "
                    f"at line {cursor + 1}: {line!r}"
                ) from exc

            if atom_id in records:
                raise ValueError(f"Duplicate PSF atom ID {atom_id} in {path}.")

            records[atom_id] = {
                "line_index": cursor,
                "atom_id": atom_id,
                "segid": tokens[1],
                "resid": tokens[2],
                "resname": tokens[3],
                "atom_name": tokens[4],
                "atom_type": tokens[5],
                "charge": charge,
                "mass": mass,
                "extra_tokens": tokens[8:],
            }
            atom_line_indices.append(cursor)

        cursor += 1

    if len(atom_line_indices) != natom:
        raise ValueError(
            f"{path} declares {natom} atoms but only {len(atom_line_indices)} "
            "atom records were parsed."
        )

    return lines, records, atom_line_indices


def validate_against_psf(summary: pd.DataFrame, psf_records: dict[int, dict[str, object]]) -> None:
    for row in summary.itertuples(index=False):
        atom_id = int(row.atom_id)
        if atom_id not in psf_records:
            raise ValueError(f"WOOPS atom ID {atom_id} is absent from the PSF.")

        record = psf_records[atom_id]
        comparisons = {
            "resid": (str(row.resid).strip(), str(record["resid"]).strip()),
            "resname": (str(row.resname).strip(), str(record["resname"]).strip()),
            "atom_name": (
                str(row.atom_name).strip(),
                str(record["atom_name"]).strip(),
            ),
        }
        for field, (woops_value, psf_value) in comparisons.items():
            if woops_value != psf_value:
                raise ValueError(
                    f"PSF mismatch for atom {atom_id}, {field}: "
                    f"WOOPS={woops_value!r}, PSF={psf_value!r}."
                )

        if not math.isclose(
            float(row.original_charge),
            float(record["charge"]),
            rel_tol=0.0,
            abs_tol=5.0e-8,
        ):
            raise ValueError(
                f"PSF charge mismatch for atom {atom_id}: "
                f"WOOPS={row.original_charge}, PSF={record['charge']}."
            )


def write_modified_psf(
    source_lines: list[str],
    psf_records: dict[int, dict[str, object]],
    values: pd.Series,
    destination: Path,
) -> int:
    lines = list(source_lines)
    replaced = 0

    for atom_id, value in values.items():
        atom_id = int(atom_id)
        if atom_id not in psf_records:
            raise KeyError(f"Cannot write atom {atom_id}; it is absent from the PSF.")

        record = psf_records[atom_id]
        line_index = int(record["line_index"])
        extra_tokens = list(record["extra_tokens"])

        new_line = (
            f"{atom_id:8d} "
            f"{str(record['segid']):<4s} "
            f"{str(record['resid']):>4s} "
            f"{str(record['resname']):<4s} "
            f"{str(record['atom_name']):<4s} "
            f"{str(record['atom_type']):<6s} "
            f"{float(value):14.10f} "
            f"{float(record['mass']):13.4f}"
        )
        if extra_tokens:
            new_line += " " + " ".join(extra_tokens)
        new_line += "\n"

        lines[line_index] = new_line
        replaced += 1

    destination.write_text("".join(lines), encoding="utf-8")
    return replaced


def count_pdb_atoms(path: Path) -> int:
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                count += 1
    return count


def make_weighted_summary(by_leg: pd.DataFrame) -> pd.DataFrame:
    by_leg = by_leg.copy()
    by_leg["weight_frames"] = by_leg["table_frames"].astype(np.int64)

    reference = (
        by_leg.sort_values(["leg", "atom_id"])
        .drop_duplicates("atom_id", keep="first")
        .set_index("atom_id")
    )
    grouped = by_leg.groupby("atom_id", sort=True)

    summary = reference[["resid", "atom_name", "resname", "original_charge"]].copy()
    summary["total_frames"] = grouped["weight_frames"].sum()
    summary["total_trajectory_ns"] = grouped["trajectory_ns"].sum()
    summary["n_legs"] = grouped["leg"].nunique()

    denominator = grouped["weight_frames"].sum().astype(float)

    for column in VALUE_COLUMNS:
        numerator = (
            (by_leg[column] * by_leg["weight_frames"])
            .groupby(by_leg["atom_id"])
            .sum()
        )
        summary[column] = numerator / denominator

        # Unweighted diagnostics across independent per-leg estimates.
        # These are not correlation-time-corrected errors.
        summary[f"{column}_leg_sd"] = grouped[column].std(ddof=1)
        summary[f"{column}_leg_sem"] = (
            summary[f"{column}_leg_sd"] / np.sqrt(summary["n_legs"])
        )

    return summary.reset_index().sort_values("atom_id").reset_index(drop=True)


def tcl_braced(value: Path | str) -> str:
    text = str(value)
    if "}" in text:
        raise ValueError(f"Cannot safely write Tcl braced value containing '}}': {text}")
    return "{" + text + "}"


def write_render_config(
    path: Path,
    root: Path,
    psf: Path,
    pdb: Path,
    output_dir: Path,
    simulation_id: str,
) -> None:
    image_dir = output_dir / "images"
    content = f"""# Auto-generated by process_woops_legs.py
set woops_root {tcl_braced(root)}
set woops_simulation_id {tcl_braced(simulation_id)}
set woops_original_psf {tcl_braced(psf)}
set woops_coordinate_pdb {tcl_braced(pdb)}
set woops_output_dir {tcl_braced(output_dir)}
set woops_image_dir {tcl_braced(image_dir)}

# 'not water' includes the 1435 protein atoms plus the integral copper atom
# in the present plastocyanin system. The mapped QuickSurf itself uses only
# protein atoms because those are the atoms present in p1p2.dat.
set woops_dipole_selection {{not water}}
set woops_base_selection {{not water}}
set woops_surface_selection {{protein}}
set woops_label_selection {{protein and resid 16 and name O}}

set woops_image_width 2000
set woops_image_height 2000

# name, modified PSF, physical minimum, physical maximum, zero-position
# within the normalized color scale.
set woops_map_specs {{
    {{p1 {tcl_braced(output_dir / 'woops_p1.psf')} -1.0 1.0 0.5}}
    {{p2 {tcl_braced(output_dir / 'woops_p2.psf')} -0.5 1.0 0.333333333333}}
    {{p21 {tcl_braced(output_dir / 'woops_p21.psf')} -0.5 0.5 0.5}}
}}
"""
    path.write_text(content, encoding="utf-8")


def main() -> int:
    args = parse_args()

    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"WOOPS root directory not found: {root}")

    simulation_id = normalize_simulation_id(
        root.name, "simulation directory name"
    )

    output_dir = (
        args.output_dir.expanduser()
        if args.output_dir is not None
        else root / "woops_processed"
    )
    if not output_dir.is_absolute():
        output_dir = root / output_dir
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    length_table_path = (
        args.length_table.expanduser()
        if args.length_table is not None
        else root.parent.parent / "trajectory_leg_lengths.dat"
    )
    if not length_table_path.is_absolute():
        length_table_path = root / length_table_path
    length_table_path = length_table_path.resolve()

    psf_path = resolve_single_file(
        explicit=args.psf,
        root=root,
        pattern="*.psf",
        description="original PSF",
    )
    pdb_path = resolve_single_file(
        explicit=args.pdb,
        root=root,
        pattern="*.pdb",
        description="coordinate PDB",
        exclude_names=("beta.pdb",),
        preferred_name="1YLB_Small_waters.pdb",
    )

    all_numbered = numbered_subdirectories(root)
    if not all_numbered:
        raise RuntimeError(f"No positive three-digit leg directories found in {root}.")

    if args.legs is None:
        requested_names = [path.name for path in all_numbered]
    else:
        requested_names = [
            normalize_three_digit(value, "leg") for value in args.legs
        ]
        requested_names = sorted(set(requested_names), key=int)

    directory_map = {path.name: path for path in all_numbered}
    missing_directories = [leg for leg in requested_names if leg not in directory_map]
    if missing_directories:
        raise FileNotFoundError(
            f"Requested leg directories do not exist: {missing_directories}"
        )

    leg_paths: list[tuple[str, Path]] = []
    missing_outputs: list[Path] = []
    for leg in requested_names:
        output_path = directory_map[leg] / args.p1p2_name
        if output_path.is_file():
            leg_paths.append((leg, output_path.resolve()))
        else:
            missing_outputs.append(output_path)

    if missing_outputs and not args.allow_missing:
        formatted = "\n".join(f"  - {path}" for path in missing_outputs)
        raise FileNotFoundError(
            "Numbered leg directories are missing p1p2.dat output:\n"
            f"{formatted}\n"
            "Use --allow-missing only when a deliberately partial average is intended."
        )

    if missing_outputs:
        for path in missing_outputs:
            print(f"WARNING: skipping missing output {path}", file=sys.stderr)

    if not leg_paths:
        raise RuntimeError("No p1p2.dat files remain to process.")

    simulation_rows = load_length_table(length_table_path, simulation_id)
    simulation_rows = simulation_rows.set_index("subdirectory", drop=False)

    frames: list[pd.DataFrame] = []
    diagnostic_rows: list[dict[str, object]] = []
    reference: pd.DataFrame | None = None

    print(f"WOOPS root:        {root}")
    print(f"Simulation:        {simulation_id}")
    print(f"PSF:               {psf_path}")
    print(f"PDB:               {pdb_path}")
    print(f"Length table:      {length_table_path}")
    print(f"Legs:              {', '.join(leg for leg, _ in leg_paths)}")
    print(f"Output directory:  {output_dir}")
    print()

    for leg, p1p2_path in tqdm(leg_paths, desc="Reading WOOPS legs", unit="leg"):
        if leg not in simulation_rows.index:
            raise ValueError(
                f"Leg {simulation_id}/{leg} is absent from {length_table_path}."
            )

        table_row = simulation_rows.loc[leg]
        if isinstance(table_row, pd.DataFrame):
            raise ValueError(
                f"Multiple length-table rows found for {simulation_id}/{leg}."
            )

        leg_df = parse_p1p2(p1p2_path, leg)
        woops_frames = int(leg_df["woops_frames"].iloc[0])
        table_frames = int(table_row["n_frames"])

        if woops_frames != table_frames:
            raise ValueError(
                f"Frame-count mismatch for {simulation_id}/{leg}: "
                f"p1p2.dat reports {woops_frames}, "
                f"trajectory_leg_lengths.dat reports {table_frames}."
            )

        if reference is None:
            reference = leg_df.copy()
        else:
            validate_leg_metadata(reference, leg_df, p1p2_path)

        leg_df["simulation_id"] = simulation_id
        leg_df["table_frames"] = table_frames
        leg_df["trajectory_ns"] = float(table_row["trajectory_ns"])
        leg_df["delta_mu_D"] = float(table_row["delta_mu_D"])
        leg_df["all_subdirs_have_dcd"] = bool(table_row["all_subdirs_have_dcd"])
        frames.append(leg_df)

        diagnostic_rows.append(
            {
                "simulation_id": simulation_id,
                "leg": leg,
                "p1p2_file": str(p1p2_path),
                "n_atoms": len(leg_df),
                "woops_frames": woops_frames,
                "table_frames": table_frames,
                "trajectory_ns": float(table_row["trajectory_ns"]),
                "delta_mu_D": float(table_row["delta_mu_D"]),
                "all_subdirs_have_dcd": bool(table_row["all_subdirs_have_dcd"]),
                "frame_count_match": True,
                "atom_metadata_match": True,
                "status": "PASS",
            }
        )

    by_leg = pd.concat(frames, ignore_index=True)
    diagnostics = pd.DataFrame.from_records(diagnostic_rows)
    summary = make_weighted_summary(by_leg)

    source_lines, psf_records, atom_line_indices = read_psf(psf_path)
    if len(atom_line_indices) != len(psf_records):
        raise RuntimeError("Internal PSF parser inconsistency.")

    pdb_atom_count = count_pdb_atoms(pdb_path)
    if pdb_atom_count != len(psf_records):
        raise ValueError(
            f"PDB/PSF atom-count mismatch: {pdb_path} has {pdb_atom_count} "
            f"ATOM/HETATM records, while {psf_path} has {len(psf_records)} atoms."
        )

    validate_against_psf(summary, psf_records)

    by_leg_columns = [
        "simulation_id",
        "leg",
        "source_file",
        "delta_mu_D",
        "trajectory_ns",
        "table_frames",
        "woops_frames",
        "all_subdirs_have_dcd",
        "atom_id",
        "resid",
        "atom_name",
        "resname",
        "original_charge",
        "p1",
        "p2",
        "p21",
    ]
    by_leg_path = output_dir / "p1p2_by_leg.dat"
    summary_path = output_dir / "p1p2_weighted_average.dat"
    diagnostics_path = output_dir / "woops_processing_diagnostics.dat"

    by_leg[by_leg_columns].to_csv(
        by_leg_path,
        sep="\t",
        index=False,
        float_format="%.10f",
        na_rep="NaN",
    )
    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
        float_format="%.10f",
        na_rep="NaN",
    )
    diagnostics.to_csv(
        diagnostics_path,
        sep="\t",
        index=False,
        float_format="%.10f",
        na_rep="NaN",
    )

    summary_indexed = summary.set_index("atom_id")
    generated_psfs: dict[str, Path] = {}
    for column in tqdm(VALUE_COLUMNS, desc="Writing mapped PSFs", unit="PSF"):
        destination = output_dir / f"woops_{column}.psf"
        replaced = write_modified_psf(
            source_lines=source_lines,
            psf_records=psf_records,
            values=summary_indexed[column],
            destination=destination,
        )
        if replaced != len(summary):
            raise RuntimeError(
                f"{destination}: replaced {replaced} atoms, expected {len(summary)}."
            )
        generated_psfs[column] = destination

    render_config_path = output_dir / "woops_render_config.tcl"
    write_render_config(
        path=render_config_path,
        root=root,
        psf=psf_path,
        pdb=pdb_path,
        output_dir=output_dir,
        simulation_id=simulation_id,
    )

    manifest = {
        "simulation_id": simulation_id,
        "root": str(root),
        "original_psf": str(psf_path),
        "coordinate_pdb": str(pdb_path),
        "trajectory_length_table": str(length_table_path),
        "legs": [leg for leg, _ in leg_paths],
        "number_of_legs": len(leg_paths),
        "number_of_mapped_atoms": int(len(summary)),
        "total_frames": int(diagnostics["table_frames"].sum()),
        "total_trajectory_ns": float(diagnostics["trajectory_ns"].sum()),
        "weighted_average_definition": "sum(n_frames_leg * p_leg) / sum(n_frames_leg)",
        "leg_sd_sem_definition": (
            "Unweighted SD and SEM across independent per-leg p estimates; "
            "not correlation-time-corrected uncertainties."
        ),
        "outputs": {
            "by_leg_table": str(by_leg_path),
            "weighted_average_table": str(summary_path),
            "diagnostics_table": str(diagnostics_path),
            "p1_psf": str(generated_psfs["p1"]),
            "p2_psf": str(generated_psfs["p2"]),
            "p21_psf": str(generated_psfs["p21"]),
            "vmd_render_config": str(render_config_path),
        },
    }
    manifest_path = output_dir / "woops_processing_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print()
    print(diagnostics.to_string(index=False))
    print()
    print(f"Processed {len(leg_paths)} legs and {len(summary)} mapped atoms.")
    print(f"Total frames: {int(diagnostics['table_frames'].sum())}")
    print(f"Total trajectory length: {diagnostics['trajectory_ns'].sum():.6f} ns")
    print()
    print("Created:")
    for path in [
        by_leg_path,
        summary_path,
        diagnostics_path,
        generated_psfs["p1"],
        generated_psfs["p2"],
        generated_psfs["p21"],
        render_config_path,
        manifest_path,
    ]:
        print(f"  {path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        raise
