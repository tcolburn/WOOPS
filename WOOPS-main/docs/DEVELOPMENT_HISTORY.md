# Development history retained during cleanup

The working archive contained several intermediate development trees. They were
used to establish lineage but intentionally excluded from the cleaned release:

- `woops_v1.1_TCmemoryfix_program`
- `woops_v1.11_TCspeedup_program`
- `woops_v1.12_TCminiSqrt_program`
- `woops_v1.13_TCwaterselect_program`
- duplicate copies of the final nanodomain source
- `woops_v3.1_1stattempt`
- pre-normalization backups such as `userDefinePlugins.f90.pre_Nfix`

The maintained source trees retained here are:

- `src/nanodomain`: final P1/P2/P21 nanodomain implementation
- `src/kirkwood`: final scalar Kirkwood time-series implementation with 1/N normalization

Compiled objects, module files, executables, editor swap files, restart/output
files, duplicated helper scripts, and project data were removed.
