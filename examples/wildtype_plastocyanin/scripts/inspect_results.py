#!/usr/bin/env python3
"""Show human-readable summaries of saved or newly generated WOOPS results."""
from pathlib import Path
import argparse
import sys
from compare_results import get_nanodomain, get_series, get_summary, SUMMARY_COLS, REGIONS
EXAMPLE=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', choices=['sample_output','runs'], default='sample_output')
    args=p.parse_args()
    r=EXAMPLE/args.results
    records=get_nanodomain(r/'nanodomain'/'p1p2.dat')
    vals=[v for _,v in records]
    n=len(vals)
    print(f'Nanodomain: {n} solute atoms, {records[0][0][0]} frames')
    for i,name in enumerate(('P1','P2','P21'), start=1):
        x=[v[i] for v in vals]
        print(f'  {name}: minimum {min(x):+.5f}, maximum {max(x):+.5f}, zeros {sum(abs(t)<1e-14 for t in x)}')
    print('  first atom:',records[0][0], 'q/P1/P2/P21 =', vals[0])
    print('Kirkwood summary (averages exclude N=0 frames):')
    for region in REGIONS:
        s=get_summary(r/'kirkwood'/f'kirkwood_{region}_summary.dat')
        d=dict(zip(SUMMARY_COLS,s))
        series=get_series(r/'kirkwood'/f'kirkwood_{region}.dat')
        print(f'  {region:12}: frames={len(series)}, <N>={d["mean_N"]:.2f}, <g_K>={d["g_k"]:.6f}, <g_perp>={d["g_perp"]:.6f}, <g_parallel>={d["g_parallel"]:.6f}, residual={d["decomposition_residual"]:+.6f}')
    return 0
if __name__=='__main__':
    sys.exit(main())
