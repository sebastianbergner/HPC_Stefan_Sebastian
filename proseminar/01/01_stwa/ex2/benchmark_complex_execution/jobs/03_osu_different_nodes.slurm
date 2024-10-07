#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name osu_different_nodes
#SBATCH --output=logs/03_osu_different_nodes.log
#SBATCH --exclusive
#SBATCH --time 10

#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=2

# load modules
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# run programs
echo "mpiexec --display-map ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency"
mpiexec --display-map ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency
echo
echo "mpiexec --display-map ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw"
mpiexec --display-map ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw