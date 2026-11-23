#!/usr/bin/env python3
"""Plot Kirkwood series from real WOOPS output (optional matplotlib dependency)."""
from pathlib import Path
import argparse
import csv
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
REGIONS=('interface','PRO16_7A','PRO16_10A')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--region',choices=REGIONS,default='PRO16_7A')
    parser.add_argument('--results',choices=['sample_output','runs'],default='sample_output')
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    path=ROOT/args.results/'kirkwood'/f'kirkwood_{args.region}.dat'
    with path.open() as handle:
        rows=list(csv.DictReader(line for line in handle if line.strip() and not line.startswith('#')))
    if not rows:
        raise ValueError(f'No rows in {path}')
    frames=[int(row['frame']) for row in rows]
    fig,ax=plt.subplots(figsize=(8.5,4.4))
    for key,label in (('g_k','g_K'),('g_perp','g_perp'),('g_parallel','g_parallel')):
        ax.plot(frames,[float(row[key]) for row in rows],marker='o',ms=4,label=label)
    ax.set_xlabel('Snapshot number (10 ps intervals)')
    ax.set_ylabel('Instantaneous Kirkwood factors')
    label = {'interface':'Full 7 Å interface', 'PRO16_7A':'PRO16, 7 Å', 'PRO16_10A':'PRO16, 10 Å'}[args.region]
    ax.set_title('WOOPS output • Wildtype plastocyanin • '+label)
    ax.legend(frameon=False)
    ax.grid(alpha=.2)
    fig.tight_layout()
    out=(args.output or ROOT/'figures'/f'{args.region}_kirkwood.png').resolve()
    out.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(out,dpi=160)
    print(f'Saved {out}')
    plt.close(fig)
if __name__=='__main__':
    main()
