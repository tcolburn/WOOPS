# WOOPS: water orientation and dipolar correlations

WOOPS analyzes how water molecules orient at a protein surface and how their dipoles correlate within a hydration region. The repository contains two Fortran programs with OpenMP parallelization.

| Analysis | Output |
| --- | --- |
| [Nanodomain maps](src/nanodomain/README.md) | Trajectory-averaged water orientation at each solute atom ($P_1$, $P_2$, $P_{21}$) |
| [Kirkwood factors](src/kirkwood/README.md) | Frame-resolved collective dipolar correlations ($g_K$, $g_{\perp}$, $g_{\parallel}$) |

## Physics and observables

Both analyses use the unit water dipole $\mathbf{u}_i$ and a local unit normal $\mathbf{n}_i$ pointing from the nearest solute heavy atom to the water oxygen. Displacements use the minimum image in an orthorhombic box. The protein interface contains all waters whose oxygens lie within 7 Å of any solute atom, including hydrogens.

### Local water orientation

Let $c_i=\mathbf{u}_i\cdot\mathbf{n}_i$, and let $\chi_i$ be the angle between the normal to the H–O–H plane and $\mathbf{n}_i\times\mathbf{u}_i$.

```math
P_1=c_i, \qquad P_2=\frac{3c_i^2-1}{2}, \qquad
P_{21}=\frac{1-c_i^2}{2}\left(2\cos^2\chi_i-1\right).
```

$P_1$ measures signed dipole alignment with the outward normal. $P_2$ distinguishes dipoles along the normal from dipoles lying in the local tangent plane. $P_{21}$ resolves the orientation of the H–O–H plane about the dipole axis and vanishes for a dipole along either direction of the normal.

Each interfacial water contributes to its nearest solute heavy atom. The program averages the assigned waters for each atom within a frame, then averages frames with equal weight. A frame with zero assigned waters contributes zero to that atom’s map. Bonded hydrogens inherit the heavy-atom values.

### Collective dipolar correlations

For a region containing $N(t)>0$ waters, define the scalar normal projection $p_i=\mathbf{u}_i\cdot\mathbf{n}_i$ and the tangential vector $\mathbf{u}_{\parallel,i}=\mathbf{u}_i-p_i\mathbf{n}_i$.

```math
\begin{aligned}
g_K(t) &= \frac{\left|\sum_i\mathbf{u}_i\right|^2}{N(t)}
       = \frac{1}{N(t)}\sum_{i,j}\mathbf{u}_i\cdot\mathbf{u}_j, \\
g_{\perp}(t) &= \frac{\left(\sum_i p_i\right)^2}{N(t)}, \\
g_{\parallel}(t) &= \frac{\left|\sum_i\mathbf{u}_{\parallel,i}\right|^2}{N(t)}.
\end{aligned}
```

The self terms give $g_K$ a baseline of 1. In the trajectory mean, positive net pair alignment raises $g_K$ above 1 and net anticorrelation lowers it. The normal factor sums signed scalar projections, so $g_{\perp}+g_{\parallel}=g_K$ for a common normal. On a curved interface the normals vary, and WOOPS reports the residual $g_K-g_{\perp}-g_{\parallel}$.

The Kirkwood executable evaluates the full 7 Å interface and two spheres of radius 7 Å and 10 Å centered on the plastocyanin PRO16 carbonyl oxygen. The site defines sphere membership, while the nearest heavy atom defines the local normal in all three regions. The trajectory summaries average over occupied frames.

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
