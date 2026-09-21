# WOOPS extensions: nanodomain maps and scalar Kirkwood analysis

This repository is a cleaned source tree containing two maintained extensions
of the WOOPS (water orientational order parameter) code originally written by
Daniel R. Martin in Dmitry Matyushov's group at Arizona State University.

> **Distribution status:** the inherited source contains an explicit 2015
> Arizona Board of Regents / ASU notice that forbids redistribution and
> derivative works without written permission. This tree is technically cleaned
> and buildable, but it should not be made public until the required permission
> or replacement license has been obtained. See `DISTRIBUTION_NOTICE.md`.

## Repository layout

```text
src/
  nanodomain/   final atom-resolved P1/P2/P21 nanodomain implementation
  kirkwood/     final scalar Kirkwood time-series implementation
tools/
  nanodomain/   P1/P2/P21 post-processing and VMD rendering helpers
  kirkwood/     scalar Kirkwood leg/time-series post-processing
examples/       placeholder for a future tiny validation/example dataset
docs/           upstream provenance and development-history notes
```

Intermediate development trees, compiled objects/modules, executables, backup
copies, editor files, and production molecular-simulation data have been
removed.

## Relationship to original WOOPS

Martin's original WOOPS workflow computes orientational order parameters for
hydration water and maps trajectory-averaged values onto a protein surface. The
2015 user guide describes PSF/PDB/DCD input, beta tagging (`beta=1` solute,
`beta=2` water), Fortran compilation, `p1p2` output, and VMD visualization. The
preserved pre-nanodomain lineage in the supplied development archive associates
the nearest water with each protein atom using a 3.0 A cutoff.

The two maintained variants here preserve the WOOPS trajectory reader and much
of its input machinery but answer different questions:

| Feature | Original WOOPS / preserved base lineage | Nanodomain variant | Scalar Kirkwood variant |
|---|---|---|---|
| Primary purpose | Atom-resolved hydration orientation map | Atom-resolved many-water P1/P2/P21 map | Collective shell dipole correlations |
| Water association | Nearest water to each protein atom in the first shell | All waters within 7 A of any protein atom; assign each to nearest protein heavy atom | Three shell selections evaluated simultaneously |
| Local direction | Protein-atom to associated water | Nearest-heavy-atom to water O | Nearest-heavy-atom to water O |
| PBC | Legacy geometry | Minimum-image local geometry | Minimum-image local geometry |
| Per-frame aggregation | One associated water per mapped atom | Average all waters assigned to each heavy atom; propagate to bonded H | Sum collective dipole/projection vectors for each region |
| Main observables | P1, P2, P21 | P1, P2, P21 | g_K, g_perp, g_parallel, residual, normal coherence |
| Main output | Trajectory-averaged atom map | Trajectory-averaged atom map | Full frame-resolved regional time series + summaries |
| Parallelization | Original serial workflow | OpenMP + periodic linked-cell search | OpenMP + periodic linked-cell search |

The nanodomain and Kirkwood programs are therefore parallel branches, not two
output modes of one executable.

## Nanodomain P1/P2/P21 branch

Source: `src/nanodomain/`

For each frame, water oxygens within 7 A of any protein atom are selected and
assigned to their nearest protein heavy atom. The water dipole is projected onto
the corresponding local surface normal, P1/P2/P21 are calculated, and all waters
assigned to a heavy atom are averaged for that frame. The heavy-atom value is
then propagated to its bonded protein hydrogens before trajectory accumulation.

See `src/nanodomain/README.md` for the exact geometry and output behavior.

## Scalar Kirkwood branch

Source: `src/kirkwood/`

This branch does not make an atom-resolved map. It evaluates the full 7 A
interface and PRO16-centered 7 A and 10 A shells during one trajectory pass. The
current corrected normalization is

```text
g_K(t)        = |sum_i u_i|^2 / N(t)
g_perp(t)     = [sum_i (u_i.n_i)]^2 / N(t)
g_parallel(t) = |sum_i [u_i - n_i(u_i.n_i)]|^2 / N(t)
```

The cleaned source intentionally contains **no `1/N^2` production
normalization**. The source comments and documentation have been updated to
match the corrected implementation.

See `src/kirkwood/README.md` for selection definitions and output columns.

## Compilation

Both maintained programs use the same simple GNU Fortran/OpenMP build pattern:

```bash
cd src/nanodomain   # or src/kirkwood
make clean
make
```

The resulting executable is `map.e`. The supplied makefiles use `gfortran -O3
-fPIC -fopenmp`. Generic 8-core Slurm launch examples are included in each
source directory.

The cleaned source trees were compilation-tested with GNU Fortran 14.2.0 during
repository preparation. This is a build check, not a substitute for the planned
independent numerical validation against a readable reference implementation.

## Python utilities

Install the small post-processing dependency set with

```bash
python -m pip install -r requirements.txt
```

The nanodomain utilities aggregate independent `p1p2.dat` legs, generate mapped
PSFs, construct difference maps, and support the VMD rendering workflow. The
Kirkwood utility combines frame-resolved leg files and writes Polars/Parquet
outputs with optional equilibration trimming.

Some nanodomain utilities originated in the plastocyanin production workflow and
therefore retain assumptions about the associated directory/metadata layout.
Production trajectories and system files are intentionally not part of this
repository.

## Water models and trajectory assumptions

The inherited reader supports three-site TIP3P (`wattype=0`) and SPCE
(`wattype=2`). TIP4P (`wattype=1`) is explicitly unsupported. The current
geometry code assumes an orthorhombic periodic box and three-site waters stored
as O,H,H in the beta=2 selection.

## Upstream reference

Daniel R. Martin and Dmitry V. Matyushov, **Dipolar Nanodomains in Protein
Hydration Shells**, *J. Phys. Chem. Lett.* **6**, 407-412 (2015).
DOI: https://doi.org/10.1021/jz5025433

See `docs/UPSTREAM_WOOPS.md` for provenance notes.

## Data and validation

Large PSF/PDB production systems, DCD trajectories, `.dat` outputs, Parquet
files, and rendered figures are deliberately excluded. Those belong in the
associated research-data archive rather than the code repository.

`examples/README.md` records the intended scope of a future tiny test fixture.
