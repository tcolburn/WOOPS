# Example-data placeholder

Production plastocyanin trajectories and analysis products are intentionally not
included in this code repository.

A future minimal example/validation package should contain only enough data to
exercise the readers and reproduce analytically controlled outputs. A useful
small fixture would include:

- one tiny PSF/topology file,
- one beta-tagged PDB,
- a one- or few-frame DCD,
- an 11-line WOOPS input file,
- expected nanodomain P1/P2/P21 values and Kirkwood values.

For Kirkwood, a particularly useful regression case is two identical unit water
dipoles in the same ensemble, for which the conventional normalization gives
`g_K = 2`. This directly detects accidental replacement of the required `1/N`
normalization by `1/N^2`.

Keep production MD data in the associated data repository (for example Zenodo),
not in Git.
