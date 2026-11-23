#!/usr/bin/env python3
"""Compare WOOPS run outputs with the bundled sample output.

This is a reproduction/compatibility check, not an independent physics
validation. See the repository's separate validation/ suite for that.
"""
from __future__ import annotations
import argparse
import csv
import math
import sys
from pathlib import Path

EXAMPLE = Path(__file__).resolve().parents[1]
REGIONS = ("interface", "PRO16_7A", "PRO16_10A")
K_COLS = ("frame", "valid", "N", "normal_coherence", "g_k", "g_perp",
          "g_parallel", "decomposition_residual")
SUMMARY_COLS = ("total_frames", "valid_frames", "zero_occupancy_frames",
                "ratio_valid", "mean_N", "valid_fraction", "mean_normal_coherence",
                "g_k", "g_perp", "g_parallel", "decomposition_residual",
                "g_perp_over_g_k", "g_parallel_over_g_k",
                "g_perp_ratio_minus_one_third", "g_parallel_ratio_minus_two_thirds")


def get_nanodomain(path: Path):
    rows = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        p = line.split()
        if len(p) != 12 or p[:3] != ['DATA', 'DI', ':']:
            raise ValueError(f'{path}:{lineno}: malformed DATA DI record')
        # 10 frames | atom ID | residue ID | atom name | residue name | q | P1 P2 P21
        meta = (int(p[3]), int(p[4]), int(p[5]), p[6], p[7])
        vals = tuple(float(s) for s in p[8:])
        if not all(math.isfinite(v) for v in vals):
            raise ValueError(f'{path}:{lineno}: nonfinite number')
        rows.append((meta, vals))
    if not rows:
        raise ValueError(f'{path}: no DATA DI records')
    if [m[1] for m, _ in rows] != list(range(1, len(rows)+1)):
        raise ValueError(f'{path}: atom IDs missing, out of order or repeated')
    return rows


def get_series(path: Path):
    with path.open(newline='') as f:
        reader = csv.DictReader(line for line in f if line.strip() and not line.lstrip().startswith('#'))
        if tuple(reader.fieldnames or ()) != K_COLS:
            raise ValueError(f'{path}: wrong columns, got {reader.fieldnames!r}')
        rows=[]
        for row in reader:
            if None in row or any(v is None for v in row.values()):
                raise ValueError(f'{path}: missing or extra CSV columns')
            vals=(int(row['frame']), int(row['valid']), int(row['N'])) + tuple(float(row[c]) for c in K_COLS[3:])
            if not all(math.isfinite(v) for v in vals):
                raise ValueError(f'{path}: nonfinite number')
            rows.append(vals)
    if not rows or [row[0] for row in rows] != list(range(1, len(rows)+1)):
        raise ValueError(f'{path}: missing, repeated or reordered frame')
    return rows


def get_summary(path: Path):
    lines=[line.split() for line in path.read_text().splitlines()
           if line.strip() and not line.lstrip().startswith('#')]
    if len(lines) != 1 or len(lines[0]) != len(SUMMARY_COLS):
        raise ValueError(f'{path}: expected one {len(SUMMARY_COLS)}-column data row')
    vals=tuple(float(v) for v in lines[0])
    if not all(math.isfinite(v) for v in vals):
        raise ValueError(f'{path}: nonfinite number')
    return vals


def close(a, b, atol, rtol):
    return math.isclose(a, b, rel_tol=rtol, abs_tol=atol)


def comparison(label, sample_output, actual, meta_indices, fields, atol, rtol, max_messages=5):
    errors=[]
    max_err=0.0
    if len(sample_output) != len(actual):
        print(f'FAIL {label}: expected {len(sample_output)} records, got {len(actual)}')
        return 1
    for i, (want, got) in enumerate(zip(sample_output, actual), 1):
        if label == 'nanodomain':
            wmeta, wvals = want
            gmeta, gvals = got
            if wmeta != gmeta:
                errors.append(f'atom row {i}: metadata {gmeta!r}, expected {wmeta!r}')
            wn, gn = wvals, gvals
        else:
            wn, gn = want, got
        for idx in meta_indices:
            if wn[idx] != gn[idx]:
                errors.append(f'row {i}, {fields[idx]}: got {gn[idx]}, expected {wn[idx]}')
        for idx in range(len(fields)):
            if idx in meta_indices:
                continue
            diff=abs(wn[idx]-gn[idx]); max_err=max(max_err,diff)
            if not close(wn[idx],gn[idx],atol,rtol):
                errors.append(f'row {i}, {fields[idx]}: got {gn[idx]:.12g}, expected {wn[idx]:.12g}, Δ={diff:.3g}')
    if errors:
        print(f'FAIL {label}: {len(errors)} discrepancies (max numerical |Δ|={max_err:.3g})')
        for x in errors[:max_messages]:
            print(f'  {x}')
        return 1
    print(f'PASS {label:30s} {len(sample_output):5d} records; max numerical |Δ|={max_err:.3g}')
    return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--atol', type=float, default=2e-9,
                        help='Absolute tolerance for output rounding differences (default 2e-9)')
    parser.add_argument('--rtol', type=float, default=1e-9,
                        help='Relative tolerance (default 1e-9)')
    parser.add_argument('--runs', type=Path, default=EXAMPLE/'runs',
                        help='Directory containing nanodomain/ and kirkwood/ actual outputs')
    args=parser.parse_args()
    if args.atol < 0 or args.rtol < 0:
        parser.error('Tolerances must be nonnegative')
    status=0
    sample_output=EXAMPLE/'sample_output'
    try:
        exp=get_nanodomain(sample_output/'nanodomain'/'p1p2.dat')
        act=get_nanodomain(args.runs/'nanodomain'/'p1p2.dat')
        status |= comparison('nanodomain',exp,act,(), ('q','P1','P2','P21'), args.atol,args.rtol)
        for region in REGIONS:
            filename=f'kirkwood_{region}.dat'
            exp=get_series(sample_output/'kirkwood'/filename)
            act=get_series(args.runs/'kirkwood'/filename)
            status |= comparison(f'Kirkwood {region}', exp, act, (0,1,2),K_COLS,args.atol,args.rtol)
            filename=f'kirkwood_{region}_summary.dat'
            exp=[get_summary(sample_output/'kirkwood'/filename)]
            act=[get_summary(args.runs/'kirkwood'/filename)]
            status |= comparison(f'Kirkwood {region} summary',exp,act,(0,1,2,3),SUMMARY_COLS,args.atol,args.rtol)
    except (ValueError, OSError, KeyError) as exc:
        print(f'FAIL cannot read/parse output: {exc}',file=sys.stderr)
        return 2
    print('OVERALL: ' + ('PASS' if status == 0 else 'FAIL'))
    return status

if __name__=='__main__':
    sys.exit(main())
