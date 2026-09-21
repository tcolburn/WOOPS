# Upstream WOOPS

The original WOOPS implementation was written by Daniel R. Martin in Dmitry
Matyushov's group at Arizona State University. The development archive supplied
for this cleanup included Martin's 2015 user guide, *Painting the Protein with
Water Orientations Surrounding Amino Acids*.

The guide describes WOOPS as a Fortran program that reads a PSF, a beta-tagged
PDB, and a DCD trajectory and computes trajectory-averaged first- and
second-rank water orientational order parameters for the first hydration layer.
The output is then mapped onto the protein surface for visualization in VMD.
The preserved pre-nanodomain source lineage uses a 3.0 A cutoff and associates
the nearest water with each protein atom.

The original scientific reference is:

Daniel R. Martin and Dmitry V. Matyushov, "Dipolar Nanodomains in Protein
Hydration Shells," *J. Phys. Chem. Lett.* **6**, 407-412 (2015).
DOI: 10.1021/jz5025433

The original user-guide PDF is not included in this cleaned tree. See
`DISTRIBUTION_NOTICE.md` before redistributing any upstream or derivative source.
