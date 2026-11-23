# Worked example: wild-type plastocyanin

The [input](input/) contains ten all-atom snapshots of wild-type plastocyanin (PDB 1YLB) in TIP3P water. The snapshots are approximately 10 ps apart and span 90 ps from first to last. The system contains 1,436 solute atoms, including copper, and 11,091 water molecules. This sample demonstrates the analysis workflow and output formats. See [MANIFEST.json](MANIFEST.json) for sampling details and [SHA256SUMS](SHA256SUMS) for file hashes.

## Run the example

From the repository root, with GNU Fortran, Make, Bash, and Python 3 installed, run

```bash
bash examples/wildtype_plastocyanin/run_example.sh
```

The script builds both programs, writes outputs and logs under [runs/](runs/), and numerically compares the outputs with [sample_output/](sample_output/). A successful comparison prints `OVERALL: PASS`. Use `--no-build` to reuse existing executables or set `OMP_NUM_THREADS` to choose a thread count. The [validation suite](../../validation/README.md) provides an independent Python comparison.

## Read the output

| Analysis | Files | What to look for |
| --- | --- | --- |
| [Nanodomain maps](../../src/nanodomain/README.md) | [`p1p2.dat`](sample_output/nanodomain/p1p2.dat) | $P_1$, $P_2$, and $P_{21}$ averaged per solute atom. Unoccupied frames contribute zero. |
| [Kirkwood factors](../../src/kirkwood/README.md) | [`kirkwood_*.dat`](sample_output/kirkwood/) | Frame-resolved factors and occupancy for each region, plus a summary over occupied frames. |

To inspect the bundled output, run the command below. Add `--results runs` to inspect a local run.

```bash
python3 examples/wildtype_plastocyanin/scripts/inspect_results.py
```

The [PRO16 7 Å plot](figures/pro16_7A_kirkwood.png) shows the three factors across the ten snapshots. To generate it with Matplotlib installed, run

```bash
python3 examples/wildtype_plastocyanin/scripts/plot_kirkwood.py --region PRO16_7A \
  --output examples/wildtype_plastocyanin/derived/pro16_7A_kirkwood.png
```

## View the nanodomain maps in VMD

From the repository root, generate the visualization PDBs and open a map in VMD.

```bash
python3 examples/wildtype_plastocyanin/scripts/make_vmd_maps.py
WOOPS_MAP=p1 vmd -e examples/wildtype_plastocyanin/scripts/view_map.tcl
```

Set `WOOPS_MAP` to `p1`, `p2`, or `p21`. The PDBs in `derived/` store map values in the temperature-factor field for coloring. Use [the beta-tagged PDB](input/wildtype.pdb) as input to WOOPS.
