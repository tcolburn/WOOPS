#!/usr/bin/env python3
"""Independent, real-trajectory validation of WOOPS nanodomain and Kirkwood.

Only depends on NumPy, SciPy and MDAnalysis, and not on WOOPS source.
Uses periodic KD-trees for hydration assignment and explicit dipole-pair
Gram matrices for the Kirkwood factors. See README.md for definitions.
"""
from __future__ import annotations

import argparse
import csv
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from MDAnalysis.coordinates.DCD import DCDReader

ROOT = Path(__file__).resolve().parent
REGIONS = ('interface', 'PRO16_7A', 'PRO16_10A')
K_FILES = {name: f'kirkwood_{name}.dat' for name in REGIONS}
K_COLS = ('frame','valid','N','normal_coherence','g_k','g_perp',
          'g_parallel','decomposition_residual')
P_COLS = ('atom_id','resid','resname','atomname','p1','p2','p21')
SUMMARY_COLS = ('total_frames','valid_frames','zero_occupancy_frames',
                'ratio_valid','mean_N','valid_fraction',
                'mean_normal_coherence','g_k','g_perp','g_parallel',
                'decomposition_residual','g_perp_over_g_k',
                'g_parallel_over_g_k','g_perp_ratio_minus_one_third',
                'g_parallel_ratio_minus_two_thirds')

@dataclass(frozen=True)
class Atom:
    aid: int
    resid: str
    resname: str
    name: str
    mass: float
    beta: int


def topology(psf: Path, pdb: Path):
    """Parse PSF atom identities/bonds and fixed-column PDB beta flags."""
    lines = psf.read_text().splitlines()
    na = next((i for i, s in enumerate(lines) if '!NATOM' in s), None)
    nb = next((i for i, s in enumerate(lines) if '!NBOND' in s), None)
    if na is None or nb is None:
        raise ValueError('PSF lacks !NATOM or !NBOND')
    num = int(lines[na].split()[0])
    records = [line.split() for line in lines[na+1:na+1+num]]
    if len(records)!=num or any(len(row)<8 for row in records):
        raise ValueError('Malformed PSF !NATOM section')
    ids=[int(row[0]) for row in records]
    if ids != list(range(1,num+1)):
        raise ValueError('PSF atom IDs must be consecutive, starting at 1')
    n_bonds=int(lines[nb].split()[0])
    tokens=[]
    for line in lines[nb+1:]:
        if '!' in line: break
        tokens.extend(line.split())
        if len(tokens)>=2*n_bonds: break
    if len(tokens)<2*n_bonds:
        raise ValueError('Incomplete PSF !NBOND section')
    bonds=[(int(tokens[2*i])-1,int(tokens[2*i+1])-1) for i in range(n_bonds)]
    pdb_lines=[s for s in pdb.read_text().splitlines() if s.startswith(('ATOM','HETATM'))]
    if len(pdb_lines)!=num:
        raise ValueError(f'PDB contains {len(pdb_lines)} atoms, PSF contains {num}')
    flags=[int(float(s[60:66])) for s in pdb_lines]
    atoms=[Atom(int(r[0]),r[2],r[3],r[4],float(r[7]),flags[i])
           for i,r in enumerate(records)]
    prot=np.array([i for i,a in enumerate(atoms) if a.beta==1],dtype=int)
    heavy=np.array([i for i in prot if atoms[i].mass>1.5],dtype=int)
    water=np.array([i for i,a in enumerate(atoms) if a.beta==2],dtype=int)
    if len(water)==0 or len(water)%3:
        raise ValueError('No water or beta=2 water atom count is not divisible by 3')
    if not np.array_equal(water,np.arange(water[0],water[0]+len(water))):
        raise ValueError('Water atoms must occupy a contiguous block')
    oxygens=water[::3]
    if any((atoms[o].name,atoms[o+1].name,atoms[o+2].name)!=('OH2','H1','H2')
           for o in oxygens):
        raise ValueError('Water topology must be ordered OH2,H1,H2')
    pro16=[i for i in prot if atoms[i].resid=='16' and atoms[i].name=='O']
    if len(pro16)!=1: raise ValueError('Expected one beta=1 PRO16 carbonyl O')
    heavy_slots={int(h):j for j,h in enumerate(heavy)}
    neighbors={int(i):[] for i in prot}
    for i,j in bonds:
        if i in neighbors: neighbors[i].append(j)
        if j in neighbors: neighbors[j].append(i)
    parents=[]
    for idx in prot:
        if idx in heavy_slots:
            parents.append(heavy_slots[int(idx)])
        else:
            candidate=[j for j in neighbors[int(idx)] if j in heavy_slots]
            if len(candidate)!=1:
                raise ValueError(f'Atom {atoms[idx].aid} lacks unique bonded heavy parent')
            parents.append(heavy_slots[candidate[0]])
    return atoms,prot,heavy,oxygens,int(pro16[0]),np.array(parents,dtype=int)


def min_image(displacement,box):
    """Orthorhombic PBC, no exact half-box ties in supplied fixture."""
    return displacement-box*np.rint(displacement/box)


def normed(v):
    length=np.linalg.norm(v,axis=-1,keepdims=True)
    if np.any(length<1e-12):raise ValueError('Degenerate vector in molecular geometry')
    return v/length


def pairwise_kirkwood(dips,normals,fast=False):
    """All-pairs sums, with optionally an algebraically equivalent fast mode."""
    n=len(dips)
    if not n:
        return dict(valid=0,N=0,normal_coherence=0.,g_k=0.,g_perp=0.,
                    g_parallel=0.,decomposition_residual=0.)
    projections=np.einsum('ij,ij->i',dips,normals)
    tangents=dips-projections[:,None]*normals
    if fast:
        # For full 500-frame local validation, avoids explicit O(N^2) matrices.
        gk=float(np.dot(np.sum(dips,axis=0),np.sum(dips,axis=0))/n)
        gp=float(np.sum(projections)**2/n)
        gt=float(np.dot(np.sum(tangents,axis=0),np.sum(tangents,axis=0))/n)
    else:
        # Mathematically independent of the Fortran's O(N) vector reduction.
        gk=float((dips@dips.T).sum()/n)
        gp=float(np.outer(projections,projections).sum()/n)
        gt=float((tangents@tangents.T).sum()/n)
    return dict(valid=1,N=n,
                normal_coherence=float(np.linalg.norm(normals.sum(axis=0))/n),
                g_k=gk,g_perp=gp,g_parallel=gt,
                decomposition_residual=gk-gp-gt)


def reference(psf,pdb,dcd,fast=False):
    atoms,prot,heavy,oxygens,pro16,parents=topology(psf,pdb)
    totals=np.zeros((len(heavy),3),dtype=float)
    regions={name:[] for name in REGIONS}
    counts=[]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',DeprecationWarning)
        with DCDReader(str(dcd)) as reader:
            if reader.n_atoms!=len(atoms):
                raise ValueError(f'DCD atom count {reader.n_atoms} != PSF {len(atoms)}')
            nframes=len(reader)
            if nframes==0:raise ValueError('Empty DCD')
            for frame,ts in enumerate(reader,1):
                angles=np.asarray(ts.dimensions[3:],dtype=float)
                if not np.allclose(angles,90.,atol=1e-3,rtol=0):
                    raise ValueError('Only orthorhombic periodic boxes supported')
                box=np.asarray(ts.dimensions[:3],dtype=float)
                if np.any(box<=0):raise ValueError('Nonpositive/missing DCD box dimensions')
                xyz=np.asarray(ts.positions,dtype=float)
                o=xyz[oxygens]
                # Two full periodic protein trees; never uses WOOPS cell lists.
                alltree=cKDTree(np.mod(xyz[prot],box),boxsize=box)
                heavy_tree=cKDTree(np.mod(xyz[heavy],box),boxsize=box)
                d_any,_=alltree.query(np.mod(o,box),k=1)
                _,nearest=heavy_tree.query(np.mod(o,box),k=1)
                nearest=np.asarray(nearest,dtype=int)
                local_normals=normed(min_image(o-xyz[heavy[nearest]],box))
                d_pro16=np.linalg.norm(min_image(o-xyz[pro16],box),axis=1)
                in_if=d_any<=7.0
                masks=dict(interface=in_if,PRO16_7A=d_pro16<=7.,PRO16_10A=d_pro16<=10.)
                relevant=in_if|(d_pro16<=10.)
                oids=oxygens[relevant]
                v1=min_image(xyz[oids+1]-xyz[oids],box)
                v2=min_image(xyz[oids+2]-xyz[oids],box)
                u=normed(v1+v2)
                n=local_normals[relevant]
                near=in_if[relevant]
                if np.any(near):
                    proj=np.clip(np.einsum('ij,ij->i',u[near],n[near]),-1,1)
                    p2=(3*proj**2-1)/2
                    water_plane=normed(np.cross(v1[near],v2[near]))
                    tangent=np.cross(n[near],u[near])
                    tangent_len=np.linalg.norm(tangent,axis=1)
                    cos2=np.zeros(len(proj),dtype=float)
                    nonzero=tangent_len>1e-12
                    cos2[nonzero]=(np.einsum('ij,ij->i',water_plane[nonzero],
                                         tangent[nonzero])/tangent_len[nonzero])**2
                    p21=(1-proj**2)*(2*cos2-1)/2
                    measures=np.column_stack((proj,p2,p21))
                    slot=nearest[in_if]
                    tally=np.bincount(slot,minlength=len(heavy))
                    frame_vals=np.zeros_like(totals)
                    for j in range(3):
                        frame_vals[:,j]=np.bincount(slot,weights=measures[:,j],
                                                    minlength=len(heavy))
                    assigned=tally>0
                    frame_vals[assigned]/=tally[assigned,None]
                    totals+=frame_vals
                counts.append(int(np.count_nonzero(in_if)))
                for name in REGIONS:
                    chosen=masks[name][relevant]
                    regions[name].append(dict(frame=frame,
                        **pairwise_kirkwood(u[chosen],n[chosen],fast=fast)))
                if frame==nframes or frame%50==0:
                    print(f'  reference {frame}/{nframes} frames; interface waters={counts[-1]}',flush=True)
    frame_mean=totals[parents]/nframes
    maps=[dict(atom_id=atoms[i].aid,resid=atoms[i].resid,resname=atoms[i].resname,
               atomname=atoms[i].name,p1=float(v[0]),p2=float(v[1]),p21=float(v[2]))
          for i,v in zip(prot,frame_mean)]
    return nframes,maps,regions,counts


def read_csv(path,key):
    with path.open(newline='') as f:
        reader=csv.DictReader(f)
        rows=list(reader)
    if not rows:raise ValueError(f'Empty baseline file: {path}')
    output={}
    for row in rows:
        k=int(row[key])
        if k in output:raise ValueError(f'Duplicate {key}={k} in {path}')
        output[k]=row
    return output


def read_nanodomain(path):
    records={}; frames=set()
    for line in path.read_text().splitlines():
        if not line.startswith('DATA DI :'):continue
        if len(line)<122:raise ValueError(f'Malformed DATA DI row: {line!r}')
        nf=int(line[10:20]); aid=int(line[20:25]);frames.add(nf)
        if aid in records:raise ValueError(f'Duplicate atom ID {aid} in {path}')
        records[aid]=dict(atom_id=aid,resid=line[28:38].strip(),
            atomname=line[38:48].strip(),resname=line[48:58].strip(),
            p1=float(line[74:90]),p2=float(line[90:106]),p21=float(line[106:122]))
    if not records:raise ValueError(f'No DATA DI rows in {path}')
    if len(frames)!=1:raise ValueError(f'Inconsistent DATA DI frame counts: {frames}')
    return frames.pop(),records


def read_kirkwood(path):
    with path.open(newline='') as f:
        rows=list(csv.DictReader(s for s in f if s.strip() and not s.lstrip().startswith('#')))
    if not rows:raise ValueError(f'No time series data in {path}')
    if not set(K_COLS).issubset(rows[0]):
        raise ValueError(f'Missing columns in Kirkwood file: {path}')
    output={}
    for row in rows:
        fid=int(row['frame'])
        if fid in output:raise ValueError(f'Duplicate frame {fid} in {path}')
        output[fid]={col:(int(row[col]) if col in ('frame','valid','N') else float(row[col]))
                     for col in K_COLS}
    return output


def read_summary(path):
    data=[s for s in path.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    if len(data)!=1:raise ValueError(f'Expected one summary data line: {path}')
    vals=data[0].split()
    if len(vals)!=len(SUMMARY_COLS):
        raise ValueError(f'{path}: expected {len(SUMMARY_COLS)} summary fields, got {len(vals)}')
    return {k:(int(v) if j<4 else float(v)) for j,(k,v) in enumerate(zip(SUMMARY_COLS,vals))}


def reference_summary(rows):
    total=len(rows)
    valid=[r for r in rows if r['valid']==1]
    nv=len(valid)
    def mean(k):return float(np.mean([r[k] for r in valid])) if nv else 0.0
    gk=mean('g_k');gp=mean('g_perp');gt=mean('g_parallel')
    ratios=abs(gk)>1e-15
    p=gp/gk if ratios else 0.
    t=gt/gk if ratios else 0.
    return dict(total_frames=total,valid_frames=nv,zero_occupancy_frames=total-nv,
        ratio_valid=int(ratios),mean_N=mean('N'),valid_fraction=nv/total,
        mean_normal_coherence=mean('normal_coherence'),g_k=gk,g_perp=gp,
        g_parallel=gt,decomposition_residual=mean('decomposition_residual'),
        g_perp_over_g_k=p,g_parallel_over_g_k=t,
        g_perp_ratio_minus_one_third=p-1/3 if ratios else 0.,
        g_parallel_ratio_minus_two_thirds=t-2/3 if ratios else 0.)


def compare_records(title,expected,observed,cols,atol,rtol,metadata=()):
    """Print diagnostic table and return True when every common key matches."""
    ok=True
    if set(expected)!=set(observed):
        print(f'FAIL {title}: missing={sorted(set(expected)-set(observed))[:8]} '
              f'extra={sorted(set(observed)-set(expected))[:8]}')
        ok=False
    common=sorted(set(expected)&set(observed))
    for col in metadata:
        mismatched=[k for k in common if str(expected[k][col])!=str(observed[k][col])]
        if mismatched:
            print(f'FAIL {title} {col}: {len(mismatched)} metadata mismatches '
                  f'(first: {mismatched[:5]})')
            ok=False
    for col in cols:
        errs=[];bad=[];rel=[]
        for k in common:
            val=float(observed[k][col]);exp=float(expected[k][col]);diff=abs(val-exp)
            errs.append(diff);rel.append(diff/max(abs(exp),1e-15))
            if not np.isclose(val,exp,atol=atol,rtol=rtol):bad.append((k,exp,val,diff))
        worst=max(errs,default=0.);worst_r=max(rel,default=0.)
        state='FAIL' if bad else 'PASS'
        print(f'{state} {title:22s} {col:26s} max_abs={worst:.4g} '
              f'max_rel={worst_r:.4g} failures={len(bad)}/{len(common)}')
        for key,e,o,delta in sorted(bad,key=lambda x:-x[3])[:4]:
            print(f'    key={key}: expected={e:.12g} observed={o:.12g} abs={delta:.4g}')
        ok &= not bad
    return ok


def compare_summary(region,expected,path,atol,rtol):
    observed=read_summary(path)
    ok=True
    for key in SUMMARY_COLS:
        x=float(observed[key]);y=float(expected[key]);error=abs(x-y)
        equal=x==y if key in SUMMARY_COLS[:4] else bool(np.isclose(x,y,atol=atol,rtol=rtol))
        if not equal:
            print(f'FAIL summary {region} {key}: Fortran={x:.12g}, Python={y:.12g} diff={error:.4g}')
            ok=False
    print(f'{"PASS" if ok else "FAIL"} {region:12s} summary (15 fields)')
    return ok


def write_reference(output,nano,regions):
    output.mkdir(parents=True,exist_ok=True)
    for name,rows,cols in [('nanodomain.csv',nano,P_COLS)]+[
            (K_FILES[reg].replace('.dat','.csv'),regions[reg],K_COLS)
            for reg in REGIONS]:
        with (output/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=cols)
            writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--psf',type=Path,default=ROOT/'input/wildtype.psf')
    p.add_argument('--pdb',type=Path,default=ROOT/'input/wildtype.pdb')
    p.add_argument('--dcd',type=Path,default=ROOT/'input/wildtype.dcd')
    p.add_argument('--nanodomain',type=Path,default=ROOT/'run/nanodomain/p1p2.dat')
    p.add_argument('--kirkwood',type=Path,default=ROOT/'run/kirkwood')
    p.add_argument('--baseline-only',action='store_true',help='Check reference against saved CSV files, do not require Fortran outputs')
    p.add_argument('--skip-baseline',action='store_true',help='Use for an alternate local trajectory, no fixture baseline check')
    p.add_argument('--write-reference',type=Path,help='Export recomputed reference CSVs to a NEW folder (does not validate Fortran)')
    p.add_argument('--fast-pairs',action='store_true',help='Use O(N) algebraic Kirkwood reduction for long custom DCDs')
    p.add_argument('--atol',type=float,default=1e-8)
    p.add_argument('--rtol',type=float,default=1e-8)
    a=p.parse_args()
    if a.atol<0 or a.rtol<0:raise ValueError('Tolerances must be nonnegative')
    if a.baseline_only and a.skip_baseline:raise ValueError('Incompatible flags: --baseline-only and --skip-baseline')
    print('WOOPS independent real-trajectory reference')
    print('PSF:',a.psf,'\nPDB:',a.pdb,'\nDCD:',a.dcd)
    nf,nano,regions,counts=reference(a.psf,a.pdb,a.dcd,fast=a.fast_pairs)
    print(f'Frames={nf}; protein atoms={len(nano)}; interface waters: '
          f'min={min(counts)}, mean={np.mean(counts):.1f}, max={max(counts)}')
    if a.write_reference:
        if a.write_reference.resolve()==(ROOT/'reference').resolve():
            raise ValueError('Refusing to overwrite committed reference/ directory')
        write_reference(a.write_reference,nano,regions)
        print('Reference exported to:',a.write_reference)
        return 0
    okay=True
    if not a.skip_baseline:
        print('\n=== Saved Python reference baseline ===')
        exp_nano={int(r['atom_id']):r for r in nano}
        base_nano=read_csv(ROOT/'reference/nanodomain.csv','atom_id')
        okay &= compare_records('baseline nanodomain',exp_nano,base_nano,('p1','p2','p21'),
                                a.atol,a.rtol,metadata=('resid','resname','atomname'))
        for reg in REGIONS:
            exp={int(r['frame']):r for r in regions[reg]}
            obs=read_csv(ROOT/'reference'/K_FILES[reg].replace('.dat','.csv'),'frame')
            okay &= compare_records('baseline '+reg,exp,obs,K_COLS[1:],a.atol,a.rtol)
    if a.baseline_only:
        print('\nOVERALL:', 'PASS' if okay else 'FAIL','(baseline-only)')
        return 0 if okay else 1
    print('\n=== Independent Python vs Fortran outputs ===')
    fortran_n,fortran_map=read_nanodomain(a.nanodomain)
    if fortran_n!=nf:
        print(f'FAIL frame count: Fortran nanodomain={fortran_n}, DCD={nf}')
        okay=False
    exp_nano={int(r['atom_id']):r for r in nano}
    okay &= compare_records('Fortran nanodomain',exp_nano,fortran_map,('p1','p2','p21'),
                            a.atol,a.rtol,metadata=('resid','resname','atomname'))
    for region in REGIONS:
        obs=read_kirkwood(a.kirkwood/K_FILES[region])
        exp={int(r['frame']):r for r in regions[region]}
        okay &= compare_records('Fortran '+region,exp,obs,K_COLS[1:],a.atol,a.rtol)
        okay &= compare_summary(region,reference_summary(regions[region]),
                         a.kirkwood/K_FILES[region].replace('.dat','_summary.dat'),a.atol,a.rtol)
    print('\nOVERALL:', 'PASS' if okay else 'FAIL')
    return 0 if okay else 1

if __name__=='__main__':
    try:sys.exit(main())
    except (FileNotFoundError,ValueError,IndexError,KeyError,UnicodeError) as exc:
        print('VALIDATION ERROR:',exc,file=sys.stderr)
        sys.exit(2)
