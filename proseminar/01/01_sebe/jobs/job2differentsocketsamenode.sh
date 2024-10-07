#!/bin/bash

#SBATCH --partition=lva
#SBATCH --job-name ps1.2
#SBATCH --output=output2differentsocket.log

#SBATCH --ntasks=2
#SBATCH --nodes=1

#SBATCH --exclusive

# load all modules
module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v # if compiling has to be done here
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55 

# normaly compiling

# executing
for i in {1..5}; do
    echo $i
    mpiexec --display-map --bind-to core --map-by socket /home/cb76/cb761014/exercise1/osu-micro-benchmarks-5.8/executables/libexec/osu-micro-benchmarks/mpi/pt2pt/osu_latency
done
for i in {1..5}; do
    echo $i
    mpiexec --display-map --bind-to core --map-by socket /home/cb76/cb761014/exercise1/osu-micro-benchmarks-5.8/executables/libexec/osu-micro-benchmarks/mpi/pt2pt/osu_bw
done

# module load hwloc/2.8.0-gcc-12.2.0-mw3u4ax
# mpiexec lstopo --of txt