#!/usr/bin/env python3
"""Create atom-resolved WOOPS difference maps from two processed simulations.

The difference convention is always

    delta = A - B

where A and B are simulation-specific WOOPS roots that already contain
``woops_processed/p1p2_weighted_average.dat`` and
``woops_processed/woops_processing_manifest.json``.

The script writes a new self-contained processing directory whose mapped PSFs
and render configuration use the same filenames and Tcl variables expected by
``render_woops_maps.tcl`` / ``render_woops_maps_surf.tcl``.

Example
-------
python make_woops_difference.py \
    --a-root /path/to/WOOPS/012 \
    --b-root /path/to/WOOPS/Wildtype \
    --output-root /path/to/WOOPS/012_minus_WT

Then render with the existing Tcl script:

cd /path/to/WOOPS/012_minus_WT
vmd -dispdev text -e ../012/render_woops_maps_surf.tcl
"""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

VALUE_COLUMNS = ("p1", "p2", "p21")
METADATA_COLUMNS = ("atom_id", "resid", "atom_name", "resname")
DEFAULT_LIMITS = {
    # Theoretical ranges of differences between two valid order parameters.
    "p1": 2.0,   # [-1, 1] - [-1, 1]
    "p2": 1.5,   # [-0.5, 1] - [-0.5, 1]
    "p21": 1.0,  # [-0.5, 0.5] - [-0.5, 0.5]
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Subtract two processed WOOPS maps atom by atom and create mapped "
            "PSFs plus a VMD render configuration. Difference = A - B."
        )
    )
    parser.add_argument(
        "--a-root",
        type=Path,
        required=True,
        help="Processed simulation root used as the minuend (A).",
    )
    parser.add_argument(
        "--b-root",
        type=Path,
        required=True,
        help="Processed simulation root used as the subtrahend (B).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help=(
            "Difference-map root. Default: a sibling directory named "
            "A_minus_B when both inputs share a parent; otherwise ./A_minus_B."
        ),
    )
    parser.add_argument(
        "--simulation-id",
        default=None,
        help="Output image/manifest identifier. Default: A_minus_B.",
    )
    parser.add_argument(
        "--reference",
        choices=("a", "b"),
        default="a",
        help=(
            "Which source provides the physical PSF/PDB used for orientation "
            "and display geometry. Default: a."
        ),
    )
    parser.add_argument(
        "--psf",
        type=Path,
        default=None,
        help="Optional explicit reference PSF, overriding the chosen manifest.",
    )
    parser.add_argument(
        "--pdb",
        type=Path,
        default=None,
        help="Optional explicit reference coordinate PDB.",
    )
    parser.add_argument(
        "--p1-limit",
        type=float,
        default=DEFAULT_LIMITS["p1"],
        help="Symmetric display limit for delta p1. Default: 2.0.",
    )
    parser.add_argument(
        "--p2-limit",
        type=float,
        default=DEFAULT_LIMITS["p2"],
        help="Symmetric display limit for delta p2. Default: 1.5.",
    )
    parser.add_argument(
        "--p21-limit",
        type=float,
        default=DEFAULT_LIMITS["p21"],
        help="Symmetric display limit for delta p21. Default: 1.0.",
    )
    parser.add_argument(
        "--copy-renderer",
        type=Path,
        default=None,
        help=(
            "Optional existing Tcl renderer to copy into OUTPUT_ROOT. "
            "The copied file keeps its original basename."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of files already present in woops_processed.",
    )
    return parser.parse_args()


def sanitize_id(text: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    cleaned = cleaned.strip("_")
    if not cleaned:
        raise ValueError(f"Could not construct a safe simulation ID from {text!r}.")
    return cleaned


def require_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{description} not found: {path}")
    return path


def load_manifest(root: Path) -> tuple[Path, dict[str, Any]]:
    root = root.expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Simulation root not found: {root}")

    manifest_path = require_file(
        root / "woops_processed" / "woops_processing_manifest.json",
        "WOOPS processing manifest",
    )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {manifest_path}: {exc}") from exc
    return manifest_path, manifest


def manifest_path_value(
    manifest: dict[str, Any],
    key: str,
    manifest_path: Path,
) -> Path:
    value = manifest.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{manifest_path} has no usable {key!r} entry.")
    return Path(value).expanduser().resolve()


def summary_path(root: Path, manifest: dict[str, Any]) -> Path:
    outputs = manifest.get("outputs", {})
    candidate = outputs.get("weighted_average_table") if isinstance(outputs, dict) else None
    if isinstance(candidate, str) and candidate.strip():
        path = Path(candidate).expanduser()
        if path.is_file():
            return path.resolve()
    return require_file(
        root / "woops_processed" / "p1p2_weighted_average.dat",
        "WOOPS weighted-average table",
    )


def load_summary(path: Path, label: str) -> pd.DataFrame:
    dataframe = pd.read_csv(
        path,
        sep="\t",
        dtype={"resid": str, "atom_name": str, "resname": str},
    )

    required = set(METADATA_COLUMNS) | {"original_charge", *VALUE_COLUMNS}
    missing = required.difference(dataframe.columns)
    if missing:
        raise ValueError(
            f"{path} ({label}) is missing required columns: {sorted(missing)}."
        )

    dataframe = dataframe.copy()
    dataframe["atom_id"] = pd.to_numeric(dataframe["atom_id"], errors="raise").astype(
        np.int64
    )
    dataframe["resid"] = dataframe["resid"].astype(str).str.strip()
    dataframe["atom_name"] = dataframe["atom_name"].astype(str).str.strip()
    dataframe["resname"] = dataframe["resname"].astype(str).str.strip()

    numeric_columns = ["original_charge", *VALUE_COLUMNS]
    for column in numeric_columns:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="raise")

    if dataframe["atom_id"].duplicated().any():
        duplicate_ids = dataframe.loc[
            dataframe["atom_id"].duplicated(keep=False), "atom_id"
        ].unique()
        raise ValueError(
            f"{path} ({label}) contains duplicate atom IDs: "
            f"{duplicate_ids[:10].tolist()}."
        )

    if not np.isfinite(dataframe[numeric_columns].to_numpy(dtype=float)).all():
        raise ValueError(f"{path} ({label}) contains NaN or infinite values.")

    bounds = {
        "p1": (-1.0, 1.0),
        "p2": (-0.5, 1.0),
        "p21": (-0.5, 0.5),
    }
    tolerance = 1.0e-8
    for column, (lower, upper) in bounds.items():
        invalid = dataframe.loc[
            (dataframe[column] < lower - tolerance)
            | (dataframe[column] > upper + tolerance),
            ["atom_id", column],
        ]
        if not invalid.empty:
            raise ValueError(
                f"{path} ({label}) contains {column} outside [{lower}, {upper}]. "
                f"Examples: {invalid.head(5).to_dict(orient='records')}"
            )

    return dataframe.sort_values("atom_id").reset_index(drop=True)


def validate_correspondence(a: pd.DataFrame, b: pd.DataFrame, a_path: Path, b_path: Path) -> None:
    if len(a) != len(b):
        raise ValueError(
            f"Mapped-atom count mismatch: {a_path} has {len(a)} rows, "
            f"whereas {b_path} has {len(b)}."
        )

    for column in METADATA_COLUMNS:
        left = a[column].astype(str).to_numpy()
        right = b[column].astype(str).to_numpy()
        mismatch = left != right
        if mismatch.any():
            index = int(np.flatnonzero(mismatch)[0])
            raise ValueError(
                f"Atom correspondence mismatch in {column!r} at sorted row {index}: "
                f"A={left[index]!r}, B={right[index]!r}; "
                f"atom IDs A/B={int(a.loc[index, 'atom_id'])}/"
                f"{int(b.loc[index, 'atom_id'])}."
            )


def optional_scalar(dataframe: pd.DataFrame, column: str) -> float | int | None:
    if column not in dataframe.columns:
        return None
    values = dataframe[column].dropna().unique()
    if len(values) != 1:
        return None
    value = values[0]
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value)
    return value.item() if hasattr(value, "item") else value


def make_difference_table(
    a: pd.DataFrame,
    b: pd.DataFrame,
    a_id: str,
    b_id: str,
) -> pd.DataFrame:
    output = a[list(METADATA_COLUMNS)].copy()
    output["source_a"] = a_id
    output["source_b"] = b_id
    output["original_charge_a"] = a["original_charge"].to_numpy(dtype=float)
    output["original_charge_b"] = b["original_charge"].to_numpy(dtype=float)
    output["delta_original_charge"] = (
        output["original_charge_a"] - output["original_charge_b"]
    )

    for column in VALUE_COLUMNS:
        output[f"{column}_a"] = a[column].to_numpy(dtype=float)
        output[f"{column}_b"] = b[column].to_numpy(dtype=float)
        output[f"delta_{column}"] = output[f"{column}_a"] - output[f"{column}_b"]

    metadata_pairs = (
        ("total_frames", "total_frames_a", "total_frames_b"),
        ("total_trajectory_ns", "total_trajectory_ns_a", "total_trajectory_ns_b"),
        ("n_legs", "n_legs_a", "n_legs_b"),
    )
    for source_column, a_column, b_column in metadata_pairs:
        a_value = optional_scalar(a, source_column)
        b_value = optional_scalar(b, source_column)
        if a_value is not None:
            output[a_column] = a_value
        if b_value is not None:
            output[b_column] = b_value

    return output


def read_psf(path: Path) -> tuple[list[str], dict[int, dict[str, Any]], list[int]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)

    natom_header_index: int | None = None
    natom: int | None = None
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

    records: dict[int, dict[str, Any]] = {}
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
                    f"Malformed numeric PSF fields in {path} at line {cursor + 1}."
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


def validate_against_psf(
    difference: pd.DataFrame,
    psf_records: dict[int, dict[str, Any]],
) -> None:
    for row in difference.itertuples(index=False):
        atom_id = int(row.atom_id)
        if atom_id not in psf_records:
            raise ValueError(f"Mapped atom ID {atom_id} is absent from the reference PSF.")

        record = psf_records[atom_id]
        comparisons = {
            "resid": (str(row.resid).strip(), str(record["resid"]).strip()),
            "resname": (str(row.resname).strip(), str(record["resname"]).strip()),
            "atom_name": (
                str(row.atom_name).strip(),
                str(record["atom_name"]).strip(),
            ),
        }
        for field, (table_value, psf_value) in comparisons.items():
            if table_value != psf_value:
                raise ValueError(
                    f"Reference PSF mismatch for atom {atom_id}, {field}: "
                    f"table={table_value!r}, PSF={psf_value!r}."
                )


def write_modified_psf(
    source_lines: list[str],
    psf_records: dict[int, dict[str, Any]],
    values: pd.Series,
    destination: Path,
) -> int:
    lines = list(source_lines)
    replaced = 0

    for atom_id, value in values.items():
        atom_id = int(atom_id)
        if atom_id not in psf_records:
            raise KeyError(f"Cannot write atom {atom_id}; absent from reference PSF.")

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


def tcl_braced(value: Path | str) -> str:
    text = str(value)
    if "}" in text:
        raise ValueError(f"Cannot safely write Tcl braced value containing '}}': {text}")
    return "{" + text + "}"


def write_render_config(
    path: Path,
    output_root: Path,
    reference_psf: Path,
    reference_pdb: Path,
    processed_dir: Path,
    simulation_id: str,
    limits: dict[str, float],
) -> None:
    image_dir = processed_dir / "images"
    content = f"""# Auto-generated by make_woops_difference.py
# Difference convention: A - B
set woops_root {tcl_braced(output_root)}
set woops_simulation_id {tcl_braced(simulation_id)}
set woops_original_psf {tcl_braced(reference_psf)}
set woops_coordinate_pdb {tcl_braced(reference_pdb)}
set woops_output_dir {tcl_braced(processed_dir)}
set woops_image_dir {tcl_braced(image_dir)}

set woops_dipole_selection {{not water}}
set woops_base_selection {{not water}}
set woops_surface_selection {{protein}}
set woops_label_selection {{protein and resid 16 and name O}}

set woops_image_width 2000
set woops_image_height 2000

# Difference maps use symmetric limits around zero, so the white point is 0.5.
set woops_map_specs {{
    {{p1 {tcl_braced(processed_dir / 'woops_p1.psf')} {-limits['p1']:.12g} {limits['p1']:.12g} 0.5}}
    {{p2 {tcl_braced(processed_dir / 'woops_p2.psf')} {-limits['p2']:.12g} {limits['p2']:.12g} 0.5}}
    {{p21 {tcl_braced(processed_dir / 'woops_p21.psf')} {-limits['p21']:.12g} {limits['p21']:.12g} 0.5}}
}}
"""
    path.write_text(content, encoding="utf-8")


def check_output_path(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists: {path}\nUse --overwrite to replace generated files."
        )


def main() -> int:
    args = parse_args()

    a_root = args.a_root.expanduser().resolve()
    b_root = args.b_root.expanduser().resolve()
    a_manifest_path, a_manifest = load_manifest(a_root)
    b_manifest_path, b_manifest = load_manifest(b_root)

    a_id = sanitize_id(str(a_manifest.get("simulation_id", a_root.name)))
    b_id = sanitize_id(str(b_manifest.get("simulation_id", b_root.name)))
    simulation_id = sanitize_id(args.simulation_id or f"{a_id}_minus_{b_id}")

    if args.output_root is None:
        parent = a_root.parent if a_root.parent == b_root.parent else Path.cwd()
        output_root = (parent / simulation_id).resolve()
    else:
        output_root = args.output_root.expanduser().resolve()

    processed_dir = output_root / "woops_processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    a_summary_path = summary_path(a_root, a_manifest)
    b_summary_path = summary_path(b_root, b_manifest)
    a_summary = load_summary(a_summary_path, "A")
    b_summary = load_summary(b_summary_path, "B")
    validate_correspondence(a_summary, b_summary, a_summary_path, b_summary_path)

    difference = make_difference_table(a_summary, b_summary, a_id, b_id)

    reference_manifest = a_manifest if args.reference == "a" else b_manifest
    reference_manifest_path = a_manifest_path if args.reference == "a" else b_manifest_path

    if args.psf is None:
        reference_psf = require_file(
            manifest_path_value(reference_manifest, "original_psf", reference_manifest_path),
            "Reference PSF",
        )
    else:
        reference_psf = require_file(args.psf, "Reference PSF")

    if args.pdb is None:
        reference_pdb = require_file(
            manifest_path_value(reference_manifest, "coordinate_pdb", reference_manifest_path),
            "Reference coordinate PDB",
        )
    else:
        reference_pdb = require_file(args.pdb, "Reference coordinate PDB")

    limits = {
        "p1": float(args.p1_limit),
        "p2": float(args.p2_limit),
        "p21": float(args.p21_limit),
    }
    for name, value in limits.items():
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"--{name}-limit must be positive and finite; received {value}.")

    difference_table_path = processed_dir / "p1p2_difference.dat"
    diagnostics_path = processed_dir / "woops_difference_diagnostics.dat"
    render_config_path = processed_dir / "woops_render_config.tcl"
    manifest_path = processed_dir / "woops_difference_manifest.json"
    generated_psf_paths = {
        column: processed_dir / f"woops_{column}.psf" for column in VALUE_COLUMNS
    }

    generated_targets = [
        difference_table_path,
        diagnostics_path,
        render_config_path,
        manifest_path,
        *generated_psf_paths.values(),
    ]
    for path in generated_targets:
        check_output_path(path, args.overwrite)

    source_lines, psf_records, atom_line_indices = read_psf(reference_psf)
    if len(atom_line_indices) != len(psf_records):
        raise RuntimeError("Internal PSF parser inconsistency.")
    pdb_atom_count = count_pdb_atoms(reference_pdb)
    if pdb_atom_count != len(psf_records):
        raise ValueError(
            f"PDB/PSF atom-count mismatch: {reference_pdb} has {pdb_atom_count} "
            f"ATOM/HETATM records, whereas {reference_psf} has {len(psf_records)} atoms."
        )
    validate_against_psf(difference, psf_records)

    difference.to_csv(
        difference_table_path,
        sep="\t",
        index=False,
        float_format="%.10f",
        na_rep="NaN",
    )

    difference_indexed = difference.set_index("atom_id")
    for column in tqdm(VALUE_COLUMNS, desc="Writing difference PSFs", unit="PSF"):
        destination = generated_psf_paths[column]
        replaced = write_modified_psf(
            source_lines=source_lines,
            psf_records=psf_records,
            values=difference_indexed[f"delta_{column}"],
            destination=destination,
        )
        if replaced != len(difference):
            raise RuntimeError(
                f"{destination}: replaced {replaced} atoms, expected {len(difference)}."
            )

    write_render_config(
        path=render_config_path,
        output_root=output_root,
        reference_psf=reference_psf,
        reference_pdb=reference_pdb,
        processed_dir=processed_dir,
        simulation_id=simulation_id,
        limits=limits,
    )

    diagnostics_rows: list[dict[str, Any]] = []
    for column in VALUE_COLUMNS:
        values = difference[f"delta_{column}"].to_numpy(dtype=float)
        diagnostics_rows.append(
            {
                "quantity": f"delta_{column}",
                "definition": f"{column}({a_id}) - {column}({b_id})",
                "n_atoms": len(values),
                "minimum": float(values.min()),
                "maximum": float(values.max()),
                "mean": float(values.mean()),
                "mean_absolute": float(np.mean(np.abs(values))),
                "rms": float(np.sqrt(np.mean(values * values))),
                "display_min": -limits[column],
                "display_max": limits[column],
                "n_below_display_min": int(np.count_nonzero(values < -limits[column])),
                "n_above_display_max": int(np.count_nonzero(values > limits[column])),
            }
        )
    diagnostics = pd.DataFrame.from_records(diagnostics_rows)
    diagnostics.to_csv(
        diagnostics_path,
        sep="\t",
        index=False,
        float_format="%.10f",
    )

    copied_renderer: Path | None = None
    if args.copy_renderer is not None:
        renderer = require_file(args.copy_renderer, "Tcl renderer")
        copied_renderer = output_root / renderer.name
        if copied_renderer.exists() and not args.overwrite:
            raise FileExistsError(
                f"Renderer destination already exists: {copied_renderer}; use --overwrite."
            )
        shutil.copy2(renderer, copied_renderer)

    manifest = {
        "simulation_id": simulation_id,
        "difference_definition": "A - B",
        "source_a": {
            "simulation_id": a_id,
            "root": str(a_root),
            "manifest": str(a_manifest_path),
            "weighted_average_table": str(a_summary_path),
        },
        "source_b": {
            "simulation_id": b_id,
            "root": str(b_root),
            "manifest": str(b_manifest_path),
            "weighted_average_table": str(b_summary_path),
        },
        "reference_geometry_source": args.reference,
        "reference_psf": str(reference_psf),
        "reference_pdb": str(reference_pdb),
        "number_of_mapped_atoms": int(len(difference)),
        "display_limits": limits,
        "outputs": {
            "difference_table": str(difference_table_path),
            "diagnostics_table": str(diagnostics_path),
            "p1_psf": str(generated_psf_paths["p1"]),
            "p2_psf": str(generated_psf_paths["p2"]),
            "p21_psf": str(generated_psf_paths["p21"]),
            "vmd_render_config": str(render_config_path),
            "copied_renderer": None if copied_renderer is None else str(copied_renderer),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print()
    print(f"Difference convention: {a_id} - {b_id}")
    print(f"Mapped atoms:          {len(difference)}")
    print(f"Reference geometry:    source {args.reference.upper()}")
    print(f"Reference PSF:         {reference_psf}")
    print(f"Reference PDB:         {reference_pdb}")
    print(f"Output root:           {output_root}")
    print()
    print(diagnostics.to_string(index=False))
    print()
    print("Created:")
    for path in generated_targets:
        print(f"  {path}")
    if copied_renderer is not None:
        print(f"  {copied_renderer}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\nERROR: {exc}", file=__import__("sys").stderr)
        raise
