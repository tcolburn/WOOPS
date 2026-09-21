#!/bin/bash
#SBATCH -J WOOPS-kirkwood
#SBATCH -o WOOPS.o%j
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -t 12:00:00

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DYNAMIC=FALSE

# Run from a directory containing the prepared 11-line WOOPS input file.
/path/to/WOOPS-extensions/src/kirkwood/map.e < WOOPS_input
