# Cleanup manifest

Source archive: `lean_WOOPS(3).zip`

The uploaded working archive contained 373 ZIP entries. The cleaned tree retains
only the two final maintained source branches, their required build files, and a
small set of useful post-processing/rendering utilities.

## Retained

- Final nanodomain P1/P2/P21 source: `src/nanodomain/`
- Final corrected scalar Kirkwood source: `src/kirkwood/`
- Nanodomain leg aggregation, difference-map, and VMD rendering helpers
- Scalar Kirkwood time-series post-processing
- Generic OpenMP Slurm launch examples
- Provenance, development-history, distribution, and example-data notes

## Removed

- Intermediate nanodomain development trees (`memoryfix`, `speedup`, `miniSqrt`, `waterselect`)
- Duplicate copies of the final nanodomain tree
- Kirkwood `v3.1` development attempt
- `userDefinePlugins.f90.pre_Nfix` backups
- All compiled `.o` and `.mod` files and `map.e` executables
- `safe/` backup directories and editor swap files
- Project-specific setup/submission scripts with hard-coded scratch paths
- Duplicate post-processing/rendering scripts
- Production PDB/PSF systems, trajectories, raw analysis outputs, logs, and images
- The upstream WOOPS user-guide PDF (provenance summarized in `docs/UPSTREAM_WOOPS.md`)

## Corrections made during cleanup

- Updated stale Kirkwood source comments from `1/N^2` to the corrected `1/N`
  normalization already present in the executable calculation.
- Replaced the stale Kirkwood README with `src/kirkwood/README.md`, documenting
  the corrected normalization and actual current outputs.
- Removed user/account-specific scratch paths from retained utility examples.
- Added `.gitignore`, Python dependency list, repository-level README, and
  distribution/provenance documentation.

## Checks performed

- `src/nanodomain` compiled successfully with GNU Fortran 14.2.0 and OpenMP.
- `src/kirkwood` compiled successfully with GNU Fortran 14.2.0 and OpenMP.
- Retained Python utilities passed `python -m py_compile`.
- No `.o`, `.mod`, `map.e`, swap, or `pre_Nfix` files remain in the packaged tree.
- No stale `N(t)^2` Kirkwood normalization remains in the maintained Kirkwood
  source/documentation.

These are build/syntax checks. Independent scientific validation of the
numerical implementation is intentionally a separate task.
