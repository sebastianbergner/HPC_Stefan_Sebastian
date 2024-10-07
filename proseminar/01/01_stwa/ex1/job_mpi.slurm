#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name mpi_test
#SBATCH --output=output_mpi.log
#SBATCH --ntasks=8
#SBATCH --nodes=8
#SBATCH --exclusive

module load openmpi/3.1.6-gcc-12.2.0-d2gmn55
mpiexec -n $SLURM_NTASKS /bin/hostname