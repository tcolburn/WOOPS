#!/usr/bin/env python3
"""Recreate the real plastocyanin validation subset from a full wildtype PSF/DCD.

This utility does not change WOOPS source, does not remove any atom, and does not
modify the checked-in reference CSVs. Always write to a new destination and
review the result before replacing the distributed validation/input/ files.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path
import warnings

import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDWriter


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--psf',type=Path,required=True,help='Original all-atom wildtype PSF')
    parser.add_argument('--dcd',type=Path,required=True,help='Last 100 ps wildtype DCD')
    parser.add_argument('--out',type=Path,required=True,help='NEW directory for generated inputs')
    parser.add_argument('--stride',type=int,default=50,help='Take one original frame every 50 frames')
    parser.add_argument('--count',type=int,default=10,help='Select ten frames by default')
    parser.add_argument('--start',type=int,default=0,help='Zero-based starting input frame')
    args=parser.parse_args()
    if args.stride<=0 or args.count<=0 or args.start<0:
        parser.error('stride/count must be positive and start must be nonnegative')
    if args.out.exists() and any(args.out.iterdir()):
        parser.error(f'Output directory is not empty; refusing overwrite: {args.out}')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',DeprecationWarning)
        system=mda.Universe(str(args.psf),str(args.dcd))
    indices=list(range(args.start,len(system.trajectory),args.stride))[:args.count]
    if len(indices)!=args.count:
        raise ValueError(f'Only {len(indices)} frames selected; requested {args.count}')
    water=np.array([a.resname=='TIP3' for a in system.atoms],dtype=bool)
    if not water.any() or water.sum()%3:
        raise ValueError('Expected three-site TIP3 water atoms')
    water_idx=np.flatnonzero(water)
    if not np.array_equal(water_idx,np.arange(water_idx[0],len(system.atoms))):
        raise ValueError('The water block must be contiguous and end with the last atom')
    wnames=system.atoms[water].names
    if any(tuple(wnames[i:i+3])!=('OH2','H1','H2') for i in range(0,len(wnames),3)):
        raise ValueError('Water atom order is not OH2,H1,H2')
    pro16=system.select_atoms('resid 16 and name O and segid P1')
    if len(pro16)!=1:raise ValueError('Expected one PRO16 P1 carbonyl oxygen')
    solute=~water
    ns=sum(solute)
    nr=len(system.atoms[solute].residues)
    args.out.mkdir(parents=True,exist_ok=True)
    # Match the legacy WOOPS pdbread fixed-column layout and beta tagging.
    original=system.trajectory[indices[0]]
    coords=np.asarray(original.positions,dtype=float).copy()
    with (args.out/'wildtype.pdb').open('w') as f:
        f.write('REMARK WOOPS WILDTYPE REAL-SYSTEM EXCERPT\n')
        for index,atom in enumerate(system.atoms):
            x,y,z=coords[index]
            tag=2.0 if water[index] else 1.0
            f.write(f'ATOM {index+1:6d} {atom.name:<4s} {atom.resname:<4s}A{int(atom.resid):4d}    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{tag:6.2f}     '
                    f'{atom.segid[:5]:<5s}{atom.name[:2]:<2s}\n')
        f.write('END\n')
    shutil.copyfile(args.psf,args.out/'wildtype.psf')
    dt=float(system.trajectory.dt)
    with DCDWriter(str(args.out/'wildtype.dcd'),n_atoms=len(system.atoms),
                   dt=dt*args.stride) as writer:
        for index in indices:
            system.trajectory[index]
            writer.write(system.atoms)
    print(f'Selected 0-based original frame indices: {indices}')
    print(f'Original dt={dt:.6g} ps; excerpt spacing={dt*args.stride:.6g} ps')
    print(f'Atoms={len(system.atoms)}; beta=1 solute={ns}; waters={water.sum()//3}; residues={nr}')
    print(f'Created input files at {args.out}')
    print('Reference CSVs are NOT updated. Run validate.py --baseline-only after replacing inputs.')


if __name__=='__main__':
    main()
