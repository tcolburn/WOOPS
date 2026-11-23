# Numerical validation

The validation suite compares the WOOPS Fortran programs with an independent Python calculation on ten snapshots of wild-type plastocyanin (PDB 1YLB) in TIP3P water. The sample retains all 34,709 atoms and uses zero-based source frames 0, 50, …, 450 from a trajectory saved every 0.2 ps. See [MANIFEST.json](MANIFEST.json) for fixture metadata and [SHA256SUMS](SHA256SUMS) for checksums.

## Run the checks

With GNU Fortran, Make, Bash, and Python 3 installed, run from the repository root.

```bash
python -m pip install -r validation/requirements.txt
bash validation/run_validation.sh
```

The runner builds both programs, writes outputs and logs under [run/](run/), and compares the results with the Python calculation. The script defaults to two OpenMP threads. Use `--no-build` to reuse existing executables or set `OMP_NUM_THREADS` to choose a thread count. To compare the Python calculation with the bundled [reference CSVs](reference/) without building Fortran, run

```bash
python validation/validate.py --baseline-only
```

## What is checked

The Python reference uses periodic SciPy `cKDTree` queries for spatial selection and evaluates explicit dipole pairs for the Kirkwood factors. The Fortran programs use linked-cell searches and collective sums. See the method descriptions for [nanodomain maps](../src/nanodomain/README.md) and [Kirkwood factors](../src/kirkwood/README.md).

| Analysis | Comparison |
| --- | --- |
| Nanodomain maps | Frame count, atom identity, and $P_1$, $P_2$, and $P_{21}$ for all 1,436 solute atoms |
| Kirkwood series | Frame index, occupancy, validity flag, three factors, normal coherence, and the residual in each region |
| Kirkwood summaries | All 15 summary fields for each region |

The default absolute and relative tolerances are both `1e-8`. Frame indices and counts are compared exactly. A successful run prints `OVERALL: PASS`. A mismatch report identifies the affected frames or atoms. This ten-snapshot test covers the spatial calculations and output conventions for the supplied orthorhombic TIP3P system.
