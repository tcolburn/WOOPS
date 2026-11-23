# WOOPS real-system validation: wildtype plastocyanin

This directory is a **drop-in numerical validation suite** for the two maintained WOOPS extensions:

- `src/nanodomain/`: atom-resolved hydration orientation maps, $P_1$, $P_2$, and $P_{21}$.
- `src/kirkwood/`: three frame-resolved collective hydration-water dipole observables and their run summaries.

All input coordinates come from the **actual wildtype plastocyanin (1YLB) simulation**. The example retains the **entire solute, copper, and solvent**; only the **number of DCD frames** has been reduced for a portable regression test. No simulated/artificial test geometry is involved.

The purpose is straightforward: compile the two maintained Fortran executables, run them on the included PSF/PDB/DCD, and compare the numerical outputs against an **independent Python implementation**. The Python script does not import Fortran code, reuse the WOOPS linked-cell implementation, or read WOOPS output when calculating its reference values.

## 1. Placement in the repository

Unzip this package so `validation/` is **directly under the repository root**, at the same level as `src/`:

```text
WOOPS-extensions/
├── README.md                    # existing repository documentation
├── src/
│   ├── nanodomain/
│   │   ├── makefile
│   │   └── ...
│   └── kirkwood/
│       ├── makefile
│       └── ...
└── validation/                  # THIS entire directory
    ├── README.md
    ├── run_validation.sh        # build, run, compare in one command
    ├── validate.py             # independent Python reference/comparison
    ├── regenerate_excerpt.py   # optional provenance/recreation utility
    ├── requirements.txt
    ├── MANIFEST.json
    ├── SHA256SUMS
    ├── input/
    │   ├── wildtype.psf         # complete, original PSF
    │   ├── wildtype.pdb         # beta-tagged, original atom order
    │   └── wildtype.dcd         # 10 full-system snapshots
    ├── reference/
    │   ├── nanodomain.csv
    │   ├── kirkwood_interface.csv
    │   ├── kirkwood_PRO16_7A.csv
    │   └── kirkwood_PRO16_10A.csv
    └── run/
        ├── nanodomain/
        │   └── WOOPS_input      # ready-to-run 11-line control file
        └── kirkwood/
            └── WOOPS_input
```

Only `validation/` is new. **Neither the Fortran source nor the existing repository files are replaced.** Run directories are intentionally included with just `WOOPS_input`; produced `.dat` files, `run.log`, and restart files are ignored by `validation/.gitignore` and should not be committed.

## 2. Quick start

Requirements:

- Python 3.10+ recommended, with NumPy, SciPy, and MDAnalysis.
- GNU Fortran with OpenMP and `make` for the WOOPS executables.
- Existing maintained `src/nanodomain` and `src/kirkwood` WOOPS branches.
- Bash for the one-command runner (`run_validation.sh`); manual commands are provided below.

From the **repository root**:

```bash
python -m pip install -r validation/requirements.txt
bash validation/run_validation.sh
```

The runner will:

1. Invoke `make` in `src/nanodomain` and `src/kirkwood`.
2. Remove **only previous generated results** from `validation/run/` to prevent a stale pass.
3. Run each WOOPS executable in its own run directory, reading `WOOPS_input` from standard input.
4. Evaluate the reference values directly from the PSF/PDB/DCD in Python.
5. Verify that the reference still matches the shipped, immutable-by-default `reference/*.csv` baselines.
6. Compare the newly produced Fortran data, atom-by-atom and frame-by-frame.
7. Check **all three Kirkwood summary files**, including mean occupancy, validity counts, ratios, and residual.
8. Print `OVERALL: PASS` or `OVERALL: FAIL` and return a corresponding shell exit code.

To reuse executables that are already compiled:

```bash
bash validation/run_validation.sh --no-build
```

The runner defaults to two OpenMP threads for this small fixture. To test another thread count, for example eight:

```bash
OMP_NUM_THREADS=8 bash validation/run_validation.sh --no-build
```

Fortran source files are not rewritten by the runner; `make` may write `.o`, `.mod`, and `map.e` build products in `src/` as usual.

### Manual run: nanodomain

```bash
make -C src/nanodomain
cd validation/run/nanodomain
../../../src/nanodomain/map.e < WOOPS_input > run.log 2>&1
cd ../../..
```

### Manual run: Kirkwood

```bash
make -C src/kirkwood
cd validation/run/kirkwood
../../../src/kirkwood/map.e < WOOPS_input > run.log 2>&1
cd ../../..
```

### Compare independently

```bash
python validation/validate.py
```

The script finds the PSF/PDB/DCD, the saved reference CSVs, and both Fortran output directories relative to `validation/validate.py`. It therefore works from *any current working directory* as long as the directory layout is preserved. To check only the Python computation and the committed baselines **without running Fortran**:

```bash
python validation/validate.py --baseline-only
```

The comparison exits with:

- `0` = passed all requested comparisons;
- `1` = numerical mismatch, missing frames/atoms, or mismatched metadata;
- `2` = bad/missing input or incompatible output format.

Example output ends with `OVERALL: PASS`. Every numerical column includes maximum absolute and relative discrepancies, plus the number of failures. When an observable fails, the script prints the first few largest discrepancies with atom IDs or frame numbers to support debugging.

## 3. Exact provenance and contents

The two files supplied as the starting point were:

- `1YLB_Small_waters.psf`: topology for the full wildtype hydrated plastocyanin system.
- `last_100ps.dcd`: **500 frames** spanning approximately 100 ps, with a 0.2 ps (200 fs) sampling interval.

For GitHub-sized validation, the original trajectory was **temporally subsampled only**, selecting the following zero-based frame indices from that 500-frame DCD:

```text
0, 50, 100, 150, 200, 250, 300, 350, 400, 450
```

Thus `validation/input/wildtype.dcd` contains **10 frames spaced approximately 10 ps apart**, covering about **90 ps between the first and last included frame**. This is **not** a continuous 10-frame, 2 ps segment; neither is it a full, undiscarded 100 ps trajectory. It is a ten-snapshot *spatial* validation sample of the last 100 ps. The DCD header records the resulting approximately 10 ps interval.

| Property | Included fixture |
|---|---:|
| Solute atoms (beta = 1), including Cu | 1,436 |
| Water molecules (TIP3P) | 11,091 |
| Water atoms (beta = 2) | 33,273 |
| Total atoms | 34,709 |
| Solute residues expected by legacy input | 100 |
| Retained frames | 10 |
| Frame separation | ~10 ps |
| Water topology order | oxygen, H1, H2 |
| DCD cell | periodic, orthorhombic |
| PRO16 backbone carbonyl O | present |

The PSF is **copied unchanged**. The beta-tagged PDB is written in the fixed-field layout expected by the original WOOPS PDB reader, using the coordinates of the first retained frame solely for PDB bookkeeping. The **actual coordinates used by both validation calculations come from the DCD**. Beta fields are 1 for solute (including Cu), 2 for TIP3P water.

The fixture is not centered or unwrapped: periodic minimum-image geometry is part of the test. No protein or water atom was removed to reduce size, so the nearest-heavy-atom search, hydration selection, and Kirkwood occupancies are tested on the complete solvent environment. See `MANIFEST.json` for fixture metadata and `SHA256SUMS` for input and reference checksums.

Verify that the files were not damaged in transit (from the `validation/` directory):

```bash
cd validation
sha256sum -c SHA256SUMS
```

On macOS, where `sha256sum` may not be installed, use `shasum -a 256 -c SHA256SUMS` instead.

## 4. What the Fortran input files contain

Both `run/nanodomain/WOOPS_input` and `run/kirkwood/WOOPS_input` use this same established 11-line input format:

```text
../../input/wildtype.psf
../../input/wildtype.pdb
../../input/wildtype.dcd
restart
100
1436
0
10
1
10
1
```

| Line | Value | Meaning |
|---|---|---|
| 1 | `../../input/wildtype.psf` | PSF relative to the run working directory |
| 2 | `../../input/wildtype.pdb` | Beta-tagged PDB |
| 3 | `../../input/wildtype.dcd` | Ten-frame DCD |
| 4 | `restart` | Restart-file basename |
| 5 | `100` | Legacy solute-residue count field |
| 6 | `1436` | Number of beta=1 solute atoms |
| 7 | `0` | TIP3P water model (`1` = unsupported TIP4P; `2` = SPC/E) |
| 8 | `10` | Number of frames to analyze |
| 9 | `1` | Start frame (one based) |
| 10 | `10` | End frame (one based) |
| 11 | `1` | WOOPS restart/leg control flag required by the current implementations |

**Do not launch either executable from the repository root** while feeding this file. The input paths are relative to `validation/run/nanodomain/` or `validation/run/kirkwood/`, respectively. `run_validation.sh` handles the working directories correctly.

## 5. Reference physics and independent calculation

For each water molecule $i$, let $\mathbf r_{Oi}$ be the oxygen position, $\mathbf r_{H1i}$ and $\mathbf r_{H2i}$ the hydrogen positions, and $\mathbf r_{ai}$ the position of its nearest *solute heavy atom* under orthorhombic minimum-image boundary conditions.

The **unit water dipole** for a neutral three-site positive-hydrogen-charge model is

$$
\mathbf u_i = \frac{\mathbf d_{H1,i}+\mathbf d_{H2,i}}
{|\mathbf d_{H1,i}+\mathbf d_{H2,i}|}, \qquad
\mathbf d_{Hk,i}=\operatorname{MIC}(\mathbf r_{Hki}-\mathbf r_{Oi}).
$$

The local outward normal is

$$
\mathbf n_i = \frac{\operatorname{MIC}(\mathbf r_{Oi}-\mathbf r_{ai})}
{|\operatorname{MIC}(\mathbf r_{Oi}-\mathbf r_{ai})|}.
$$

For this validation, water-selection membership is based on **distance to any beta=1 solute atom, including hydrogens**; selected waters are assigned to the **nearest solute heavy atom**, which may be different from their nearest solute atom. The two tests are deliberately kept separate. The Python implementation uses periodic `scipy.spatial.cKDTree` queries, whereas the production Fortran uses its own OpenMP linked-cell search.

### Nanodomain $P_1,P_2,P_{21}$

Define

$$
c_i = \mathbf u_i\cdot\mathbf n_i, \qquad
P_1(i)=c_i, \qquad P_2(i)=\frac{3c_i^2-1}{2}.
$$

For $P_{21}$, let $\mathbf w_i$ be the unit normal to the water H–O–H plane, and let $\chi_i$ be its angle to the normalized vector $\mathbf n_i\times\mathbf u_i$ when the latter is nonzero. Then the implemented orientation parameter is

$$
P_{21}(i)=\frac{1-c_i^2}{2}\left(2\cos^2\chi_i-1\right).
$$

A zero-length $\mathbf n_i\times\mathbf u_i$ is treated by a zero contribution after applying the prefactor $1-c_i^2$. The normal to the water plane may be taken with either cross-product orientation because $P_{21}$ uses $\cos^2\chi_i$.

The map is **not** a single mean over all water observations. In each frame, for every solute heavy atom $a$, the selected waters assigned to it are averaged:

$$
\overline{P}_\ell(a,t)=
\begin{cases}
N_a(t)^{-1}\sum_{i\to a}P_\ell(i,t),&N_a(t)>0,\\
0,&N_a(t)=0.
\end{cases}
$$

The per-frame mean is then accumulated using **all frames**, including frames with zero assigned waters:

$$
P_\ell(a)=\frac{1}{N_f}\sum_{t=1}^{N_f}\overline{P}_\ell(a,t),
\qquad \ell\in\{1,2,21\}.
$$

The heavy-atom value is copied to each protein hydrogen bonded to that heavy atom, preserving the atom-resolved output convention. The test checks the entire 1,436-solute-atom `DATA DI :` table in the Fortran `p1p2.dat`, including atom identity and the declared number of frames.

### Kirkwood observables and three water regions

The current Kirkwood implementation computes **three ensembles in one trajectory pass**:

1. `interface`: water O within 7 Å of any beta=1 protein/solute atom (including H).
2. `PRO16_7A`: water O within 7 Å of the backbone carbonyl oxygen of residue PRO16.
3. `PRO16_10A`: water O within 10 Å of that same PRO16 carbonyl oxygen.

For every selected water in **all three** ensembles, the local direction $\mathbf n_i$ is still taken relative to its nearest protein **heavy atom**. The PRO16 carbonyl is a region-selection center, *not* a shared projection axis.

For the instantaneous water count $N(t)$, let $p_i=\mathbf u_i\cdot\mathbf n_i$ and $\mathbf u_{\parallel,i}=\mathbf u_i-p_i\mathbf n_i$. Then

$$
\begin{aligned}
g_K(t) &= \frac{1}{N(t)}\sum_{i,j}\mathbf u_i\cdot\mathbf u_j,\\[3pt]
g_{\perp}(t) &= \frac{1}{N(t)}\sum_{i,j}p_i p_j,\\[3pt]
g_{\parallel}(t) &= \frac{1}{N(t)}\sum_{i,j}\mathbf u_{\parallel,i}\cdot\mathbf u_{\parallel,j}.
\end{aligned}
$$

The production code evaluates the equivalent $O(N)$ collective-vector reductions, but Python's **default validation computes the explicit $i,j$ Gram matrices**, avoiding the production reduction as a numerical implementation. An optional `--fast-pairs` flag uses the equivalent $O(N)$ formulas for longer local trajectories.

The other tested quantities are

$$
\mathrm{normal\_coherence}(t)=\frac{|\sum_i\mathbf n_i|}{N(t)},
\qquad
R(t)=g_K(t)-g_\perp(t)-g_\parallel(t).
$$

Because the individual normals vary along the protein surface, **$R(t)$ is not generally zero**. It would be incorrect for the test to impose a common-normal decomposition identity on these curved hydration-shell regions. The residual itself is compared as a numerical observable.

The comparison covers the `frame`, `valid`, `N`, `normal_coherence`, `g_k`, `g_perp`, `g_parallel`, and `decomposition_residual` columns of each region's `.dat` time series. `valid` and `N` are integer-exact comparisons. The reference also reconstructs the **15-field Kirkwood summary**, averaging only frames where $N(t)>0$, and compares all its integer and floating-point fields to the corresponding `kirkwood_*_summary.dat` output.

## 6. Recorded validation result

The distributed reference CSVs were computed independently from the same 10-frame real-system DCD, and both current WOOPS branches were run on the fixture. A comparison completed successfully with NumPy/SciPy/MDAnalysis, GNU Fortran, and two OpenMP threads.

| Output | Comparison scope | Largest observed absolute difference |
|---|---|---:|
| Nanodomain | $P_1,P_2,P_{21}$, all 1,436 solute atoms | $5.0\times10^{-11}$ |
| Kirkwood interface | Frame-resolved components/residual | $7.5\times10^{-14}$ |
| Kirkwood PRO16 7 Å | Frame-resolved components/residual | $5.3\times10^{-15}$ |
| Kirkwood PRO16 10 Å | Frame-resolved components/residual | $4.0\times10^{-15}$ |
| Kirkwood summaries | All 15 fields, all three regions | PASS |

The nanodomain file prints ten decimal places; numerical differences at roughly $10^{-11}$ are consistent with formatting. The regression defaults to `atol=rtol=1e-8`, comfortably tighter than typical scientific-analysis tolerances but looser than the formatted-output rounding scale.

To adjust these deliberately:

```bash
python validation/validate.py --atol 1e-9 --rtol 1e-9
```

The saved `reference/*.csv` files are **not overwritten by ordinary validation runs**. If you need to export newly computed reference tables, specify a *different* directory:

```bash
python validation/validate.py --write-reference /tmp/woops_reference_recomputed
```

## 7. Recreate the 10-frame fixture (optional)

If you still have the original 500-frame/100 ps `last_100ps.dcd` and `1YLB_Small_waters.psf`, you can independently reproduce the sampling and beta-tagged PDB, **without changing the committed input or baseline files**:

```bash
python validation/regenerate_excerpt.py \
  --psf /path/to/1YLB_Small_waters.psf \
  --dcd /path/to/last_100ps.dcd \
  --out /tmp/woops_new_validation_input \
  --start 0 --stride 50 --count 10
```

The utility refuses to overwrite a nonempty output directory. Check its output, then only deliberately replace the committed `validation/input/` files if intended. Any changed coordinates or sampling will cause the baseline check to fail until you separately inspect and intentionally update reference data.

For a full **500-frame local validation**, retain the all-atom PSF/PDB, point the WOOPS input files to the longer DCD, and update input lines 8 and 10 from `10` to `500`. Then compare with the corresponding custom `--dcd` file and skip the ten-frame baseline check:

```bash
python validation/validate.py \
  --dcd /path/to/last_100ps.dcd \
  --nanodomain /path/to/full_nanodomain/p1p2.dat \
  --kirkwood /path/to/full_kirkwood \
  --skip-baseline --fast-pairs
```

`--fast-pairs` changes *only the independent Python method of evaluating the same double sums* to an algebraically equivalent collective-vector reduction; it is useful for a 500-frame local comparison. It is **off by default** for the public 10-frame fixture.

## 8. What a pass does—and does not—establish

A pass checks the current Fortran numerical implementation against Python on a **real, fully hydrated protein geometry**, including the protein/water selection, periodic minimum-image distances, normal vectors, water dipoles, projection geometry, per-frame averaging, atom map inheritance, frame-resolved Kirkwood quantities, and the aggregate summary format.

It **does not** by itself validate alternative water models, nonorthorhombic boxes, unusual topology orderings, all possible edge cases, the correctness of simulation force-field physics, physical relaxation rates, 100 ps temporal averages, multiple-leg post-processing, or restart/recovery behavior. Ten widely spaced frames are appropriate for a portable spatial-algorithm regression; they are **not** appropriate for estimating time-correlation functions or relaxation times.

A mismatch should be investigated rather than eliminated by raising tolerances: inspect `run/*/run.log`, the first reported offending atom IDs or frames, the water/solute beta-tagging, the 11-line input, the DCD box, and whether the Fortran executable was compiled from the expected source branch. If the fixture itself was intentionally changed, compare its new independent reference before updating the committed baseline.

## 9. GitHub recommendations

Commit **everything under `validation/` except generated run products** (`p1p2.dat`, `kirkwood_*.dat`, logs, restart files, and caches). The 10-frame DCD is approximately 4 MiB, well under GitHub's normal per-file limit. The local `validation/.gitignore` explicitly allows the three small fixture inputs even if the repository root ignores all `.psf`, `.pdb`, and `.dcd` files. The original ~199 MiB 500-frame DCD should **not** be committed here: retain it in your local data archive.

This directory is deliberately standalone. It does not require modifying the repository README, Python post-processing utilities, or source branches. To make it visible to newcomers, you may add a one-line link from the root README to `validation/README.md` after merging.
