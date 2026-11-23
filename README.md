# WOOPS: water orientation and dipolar correlations

WOOPS analyzes hydration water in molecular dynamics trajectories. The repository contains two Fortran programs with OpenMP parallelization.

| Analysis | Output |
| --- | --- |
| [Nanodomain maps](src/nanodomain/README.md) | Atom-resolved averages of the water orientational order parameters $P_1$, $P_2$, and $P_{21}$ |
| [Kirkwood factors](src/kirkwood/README.md) | Frame-resolved dipole correlations for the 7 Å protein interface and the 7 Å and 10 Å regions centered on the plastocyanin PRO16 carbonyl oxygen |

## Run the plastocyanin example

With Bash, GNU Fortran, Make, and Python 3 installed, run from the repository root.

```bash
bash examples/wildtype_plastocyanin/run_example.sh
```

The script builds both executables, analyzes ten snapshots of wild-type plastocyanin in TIP3P water, and compares the results with the bundled sample output. The [example guide](examples/wildtype_plastocyanin/README.md) describes the output and visualization. The [validation suite](validation/README.md) compares both programs with an independent Python calculation.

For your own trajectories, see the [input format](src/nanodomain/README.md#input) and the [Kirkwood selection requirements](src/kirkwood/README.md#regions). Install [post-processing dependencies](requirements.txt) with `python -m pip install -r requirements.txt` to use the [nanodomain](tools/nanodomain/) and [Kirkwood](tools/kirkwood/) utilities.

## Citation and distribution

WOOPS was developed by Daniel R. Martin in Dmitry V. Matyushov’s group at Arizona State University. See [Martin and Matyushov, *J. Phys. Chem. Lett.* **6**, 407–412 (2015)](https://doi.org/10.1021/jz5025433).

See the [source provenance](docs/UPSTREAM_WOOPS.md) and the [distribution notice](DISTRIBUTION_NOTICE.md) for attribution and redistribution details.
