# WOOPS nanodomain P1/P2/P21 variant

This is the maintained atom-resolved WOOPS extension used to map local hydration
orientation around a protein surface.

## What changed from the original WOOPS workflow

The preserved upstream lineage associated a nearest water with each protein atom
inside a short first-shell cutoff. This variant instead treats the local
hydration environment as a many-water nanodomain on each frame:

1. Select every water oxygen within 7 A of any beta=1 protein atom, including H.
2. Assign each selected water uniquely to its nearest beta=1 protein **heavy** atom.
3. Define the local normal from that heavy atom to the water oxygen using minimum-image PBC.
4. Form the normalized three-site water dipole and evaluate

   `P1 = u.n`

   `P2 = (3 (u.n)^2 - 1)/2`

   `P21 = ((1-(u.n)^2) (2 cos^2(chi_w)-1))/2`.

5. Average all waters assigned to the same heavy atom in that frame.
6. Propagate the heavy-atom frame value to its bonded protein hydrogens.
7. Accumulate one atom-level value per frame and write trajectory averages.

Frames in which a heavy atom receives no assigned water contribute zero to the
ordinary trajectory map average, preserving the map-style normalization used by
WOOPS.

The implementation uses a periodic linked-cell search and OpenMP thread-local
accumulation. The water selection is still defined geometrically, not by the
cell partition itself.

## Build

GNU Fortran with OpenMP:

```bash
make clean
make
```

The executable is `map.e`.

## Input

The program keeps the established 11-line WOOPS input format:

1. PSF filename
2. beta-tagged PDB filename
3. DCD filename
4. restart basename
5. residue-count field
6. number of beta=1 protein atoms
7. water type (`0=TIP3P`, `1=TIP4P` unsupported, `2=SPCE`)
8. number of frames to analyze
9. start frame
10. end frame
11. restart flag

## Output

The principal output is `p1p2.dat`. Each `DATA DI` row contains the analyzed
frame count, atom metadata, original partial charge, and the trajectory-averaged
P1, P2, and P21 values.

## Utilities

See `../../tools/nanodomain/` for leg aggregation, difference-map generation,
and the VMD surface renderer used in the plastocyanin workflow.
