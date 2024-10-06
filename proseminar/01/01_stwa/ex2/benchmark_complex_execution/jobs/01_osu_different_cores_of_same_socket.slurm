#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name osu_different_cores_of_same_socket
#SBATCH --output=logs/01_osu_different_cores_of_same_socket.log
#SBATCH --exclusive
#SBATCH --time 10

#SBATCH --ntasks=2
#SBATCH --nodes=1

# load modules
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

# run programs
echo "mpiexec --display-map -bind-to core --map-by core ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency"
mpiexec --display-map --bind-to core --map-by core ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_latency
echo
echo "mpiexec --display-map -bind-to core --map-by core ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw"
mpiexec --display-map --bind-to core --map-by core ../osu-micro-benchmarks-5.8/mpi/pt2pt/osu_bw