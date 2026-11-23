# Kirkwood factors

## Regions

The executable evaluates three regions in one trajectory pass.

| Output suffix | Water oxygen selection |
| --- | --- |
| `interface` | Within 7 Å of any solute atom tagged `1`, including hydrogens |
| `PRO16_7A` | Within 7 Å of the PRO16 carbonyl oxygen |
| `PRO16_10A` | Within 10 Å of the PRO16 carbonyl oxygen |

The site lookup requires exactly one solute atom with PSF residue ID `16` and atom name `O`. For each selected water, the local unit normal $\mathbf{n}_i$ points from its nearest solute heavy atom to the water oxygen under the minimum-image convention. The carbonyl oxygen sets membership in the two PRO16 regions.

## Observables

For unit dipoles $\mathbf{u}_i$ and instantaneous occupancy $N(t)>0$, let $p_i=\mathbf{u}_i\cdot\mathbf{n}_i$ and $\mathbf{u}_{\parallel,i}=\mathbf{u}_i-p_i\mathbf{n}_i$.

```math
\begin{aligned}
g_K(t)&=\frac{\left|\sum_i\mathbf{u}_i\right|^2}{N(t)},\\
g_{\perp}(t)&=\frac{\left(\sum_i p_i\right)^2}{N(t)},\\
g_{\parallel}(t)&=\frac{\left|\sum_i\mathbf{u}_{\parallel,i}\right|^2}{N(t)}.
\end{aligned}
```

These reductions evaluate the corresponding dipole pair sums, including self terms, in $O(N)$ operations. The output also includes normal coherence $\left|\sum_i\mathbf{n}_i\right|/N(t)$ and the residual $g_K-g_{\perp}-g_{\parallel}$. The residual vanishes for a common normal and can be nonzero when the normals vary across a region.

## Run and output

Use the [WOOPS input format](../nanodomain/README.md#input) with line 11 set to `1`. For the bundled input, run these commands from the repository root.

```bash
make -C src/kirkwood
(cd examples/wildtype_plastocyanin/runs/kirkwood &&
  OMP_NUM_THREADS=2 ../../../../src/kirkwood/map.e < WOOPS_input)
```

Each region produces a CSV-formatted `kirkwood_<region>.dat` with the following header and a `kirkwood_<region>_summary.dat` with run statistics.

```csv
frame,valid,N,normal_coherence,g_k,g_perp,g_parallel,decomposition_residual
```

The `frame` column counts analyzed frames starting at 1. Empty regions have `valid=0` for that frame. Summary means are computed over frames with `valid=1`.

Analyze trajectory legs separately. The [Kirkwood post-processor](../../tools/kirkwood/process_kirkwood_timeseries.py) combines numbered legs into Parquet time series and uses `trajectory_leg_lengths.dat` to reconstruct elapsed time.
