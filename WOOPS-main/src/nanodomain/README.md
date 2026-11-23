# Nanodomain maps

WOOPS maps the orientation of hydration water onto solute atoms. In each frame, it selects water oxygens within 7 Å of any solute atom, including hydrogens, and assigns each water to its nearest solute heavy atom. The local unit normal $\mathbf{n}$ points from that heavy atom to the water oxygen. Displacements use the orthorhombic minimum image.

For the unit water dipole $\mathbf{u}$, set $c=\mathbf{u}\cdot\mathbf{n}$. Let $\chi$ be the angle between the normal to the H–O–H plane and $\mathbf{n}\times\mathbf{u}$.

```math
P_1=c, \qquad P_2=\frac{3c^2-1}{2}, \qquad
P_{21}=\frac{1-c^2}{2}\left(2\cos^2\chi-1\right).
```

For parallel or antiparallel $\mathbf{u}$ and $\mathbf{n}$, $P_{21}=0$. For each heavy atom, the program averages over assigned waters within a frame and then gives each analyzed frame equal weight. An unoccupied frame contributes zero. Bonded solute hydrogens inherit their heavy atom’s values.

## Input

Use a PSF, a beta-tagged PDB in the [supplied fixed-column format](../../examples/wildtype_plastocyanin/input/wildtype.pdb), and a DCD with an orthorhombic cell record in every frame. The three files must use the same atom order. The PDB reader skips the first line and uses the temperature factor to select solute (`1`), water (`2`), and excluded atoms (`0`). Each selected solute hydrogen must have exactly one selected heavy-atom parent in the PSF. The water selection must be ordered O, H, H by molecule.

Both executables read the same positional control file from standard input. The [example control file](../../examples/wildtype_plastocyanin/runs/nanodomain/WOOPS_input) shows the format.

| Line | Value |
| --- | --- |
| 1–3 | PSF, PDB, and DCD paths, one per line |
| 4 | Restart file basename. Use `restart` for the example. |
| 5 | Number of solute residues, retained as an input field |
| 6 | Number of solute atoms tagged `1` |
| 7 | Water model. Use `0` for TIP3P (`TIP3` in the PSF) or `2` for SPC/E (`SPCE` in the PSF). |
| 8 | Number of frames to analyze |
| 9–10 | First and last DCD frames, one-based and inclusive |
| 11 | Use `1` for a standalone run. Kirkwood requires `1`. |

Line 8 must equal the last frame minus the first frame plus one. Relative paths are resolved from the run directory. The reader accepts up to 60 characters per path or restart filename.

## Run and output

Build with GNU Fortran and OpenMP, then run from the directory containing `WOOPS_input`. For the bundled input, run these commands from the repository root.

```bash
make -C src/nanodomain
(cd examples/wildtype_plastocyanin/runs/nanodomain &&
  OMP_NUM_THREADS=2 ../../../../src/nanodomain/map.e < WOOPS_input)
```

The program writes `p1p2.dat` in the run directory. Each `DATA DI :` row contains the frame count for a standalone run, the PSF atom ID, residue ID, atom name, residue name, partial charge, and the three dimensionless order parameters.

See the [nanodomain tools](../../tools/nanodomain/) for trajectory-leg aggregation, difference maps, and VMD rendering.
