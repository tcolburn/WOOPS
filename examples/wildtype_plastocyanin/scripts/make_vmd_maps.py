#!/usr/bin/env python3
"""Create three visualization-only PDBs with WOOPS map values in B-factor.

These derived PDBs are *not* suitable as WOOPS inputs because beta=1/2
solute/water selection flags are overwritten by order-parameter values.
"""
from __future__ import annotations
from pathlib import Path
import argparse
from compare_results import get_nanodomain
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', choices=['sample_output','runs'], default='sample_output')
    args=parser.parse_args()
    data=get_nanodomain(ROOT/args.results/'nanodomain'/'p1p2.dat')
    values={m[1]:v for m,v in data}
    source=ROOT/'input/wildtype.pdb'
    pdb_lines=source.read_text().splitlines(keepends=True)
    atom_lines=[line for line in pdb_lines if line.startswith(('ATOM  ','HETATM'))]
    count=len(atom_lines)
    if count != 34709 or len(values) != 1436:
        raise ValueError(
            f'Expected 34,709 PDB atoms and 1,436 solute records, '
            f'found {count:,} PDB atoms and {len(values):,} solute records'
        )
    out=ROOT/'derived'
    out.mkdir(exist_ok=True)
    for component,label in [(1,'p1'),(2,'p2'),(3,'p21')]:
        output=[]
        for line in pdb_lines:
            if not line.startswith(('ATOM  ','HETATM')):
                output.append(line)
                continue
            aid=int(line[6:11]); val=values[aid][component] if aid in values else 0.0
            output.append(f'{line[:60]}{val:6.2f}{line[66:]}')
        path=out/f'{label}_map.pdb'
        path.write_text('REMARK BETA COLUMN HOLDS WOOPS '+label.upper()+' VALUES; NOT A VALID WOOPS INPUT\n'+''.join(output))
        print(f'Wrote {path.relative_to(ROOT)}')
    print('Use input/wildtype.pdb for analysis and the derived PDBs for visualization.')
if __name__=='__main__':main()
