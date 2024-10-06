#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name osu_different_sockets_of_same_node
#SBATCH --output=logs/02_osu_different_sockets_of_same_node_%j.log
# #SBATCH --open-mode=append
#SBATCH --exclusive
#SBATCH --time 10

#SBATCH --ntasks=2
#SBATCH --nodes=1

# load modules
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# run programs
echo "mpiexec --display-map --bind-to core --map-by socket ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency"
mpiexec --display-map --bind-to core --map-by socket ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency
echo
echo "mpiexec --display-map --bind-to core --map-by socket ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw"
mpiexec --display-map --bind-to core --map-by socket ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw