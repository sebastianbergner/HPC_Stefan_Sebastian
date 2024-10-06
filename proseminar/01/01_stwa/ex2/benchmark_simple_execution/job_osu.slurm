#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name osu_test
#SBATCH --output=output_osu.log
#SBATCH --exclusive
#SBATCH --time 10

#SBATCH --ntasks=2
#SBATCH --nodes=1

# load modules
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# run programs
echo "mpiexec ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency"
mpiexec ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency
echo
echo "mpiexec ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw"
mpiexec ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw