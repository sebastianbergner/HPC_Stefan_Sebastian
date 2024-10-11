#!/bin/bash

#SBATCH --partition=lva
#SBATCH --job-name ps2.1
#SBATCH --output=output_ex02.log

#SBATCH --ntasks=96
#SBATCH --nodes=8

#SBATCH --exclusive

# load all modules
module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v # if compiling has to be done here
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55 

# normaly compiling
mpicc monte_carlo_serial.c -O2 -o monte_carlo_serial
mpicc monte_carlo_parallel.c -O2 -o monte_carlo_parallel

# serial
for N in 1 10 100 1000 10000 100000 1000000 10000000 100000000 1000000000; do
    for i in {1..5}; do
        echo $i
        ./monte_carlo_serial $N
    done

    # parallel
    for i in {1..5}; do
        echo $i
        mpiexec monte_carlo_parallel $N # --display-map has been omited as it spams the log
    done
done
