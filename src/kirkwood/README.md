# WOOPS scalar Kirkwood variant

This is the maintained scalar Kirkwood extension of the streaming/OpenMP WOOPS
reader. It intentionally does **not** calculate atom-resolved P1/P2/P21 maps.
Instead it evaluates collective dipolar correlations for three water ensembles
in one DCD pass.

## Water ensembles

For every analyzed frame:

1. `interface_7A`: water O within 7 A of any beta=1 protein atom, including H.
2. `PRO16_7A`: water O within 7 A of the PRO16 carbonyl O.
3. `PRO16_10A`: water O within 10 A of the PRO16 carbonyl O.

The PRO16 cutoffs determine membership only. For every selected water, the local
normal is still the minimum-image vector from the nearest protein heavy atom to
the water oxygen.

## Implemented observables

For selected water `i`, let `u_i` be the unit water dipole, `n_i` its local
normal, `p_i = u_i.n_i`, and

`u_i,parallel = u_i - n_i (u_i.n_i)`.

For instantaneous occupancy `N(t)`, the production code uses the conventional
Kirkwood normalization

```text
g_K(t)        = |sum_i u_i|^2 / N(t)
g_perp(t)     = [sum_i p_i]^2 / N(t)
g_parallel(t) = |sum_i u_i,parallel|^2 / N(t)
```

These reduced expressions are algebraically equivalent to their explicit
`i,j` double sums with the same `1/N(t)` normalization, while requiring O(N)
work. Frames with `N=0` are retained in the compact time series with `valid=0`
but excluded from trajectory averages.

Because the local normal generally varies from water to water, the code also
records

`decomposition_residual = g_K - g_perp - g_parallel`

and

`normal_coherence = |sum_i n_i| / N`.

The common-normal identity `g_K = g_perp + g_parallel` is therefore a diagnostic,
not a constraint on the curved full interface.

## Internal checks

Before trajectory analysis, the executable checks the reduced O(N) formulas
against explicit `i,j` sums and checks the common-normal/isotropic 1/3 : 2/3
limits. These checks use the same `1/N` normalization as the production code.

## Build

```bash
make clean
make
```

The executable is `map.e`.

## Input

The established 11-line WOOPS input format is retained. For this dedicated
Kirkwood executable, input line 11 must be `1`; each DCD leg is analyzed
independently and combined in post-processing.

## Output

Compact frame-resolved CSV-formatted files:

- `kirkwood_interface.dat`
- `kirkwood_PRO16_7A.dat`
- `kirkwood_PRO16_10A.dat`

Columns:

`frame, valid, N, normal_coherence, g_k, g_perp, g_parallel, decomposition_residual`

A summary file is written for each region. Optional diagnostic reduction sums
can be enabled in `modules.f90` with `writeDiagnosticSums=.TRUE.`.

The full time series is enabled by default with `writeTimeSeries=.TRUE.` and
`timeSeriesStride=1`.

## Post-processing

`../../tools/kirkwood/process_kirkwood_timeseries.py` combines numbered DCD-leg
outputs, reconstructs cumulative time from `trajectory_leg_lengths.dat`, writes
combined Parquet time series, and supports an initial-time discard for summary
statistics without modifying the raw Fortran output.
