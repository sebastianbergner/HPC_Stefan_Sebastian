#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps2.2
#SBATCH --output=output/output.log

#SBATCH --ntasks=2
# #SBATCH --ntasks=96
# #SBATCH --nodes=8

#SBATCH --exclusive

# load all modules
module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v # if compiling has to be done here
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55 

make --no-print-directory

REP=1
# PROBLEMS="768 1536 3072 6144"
PROBLEMS="1536"

for i in $(seq 1 $REP); do
    for N in $PROBLEMS; do
        echo ----------------- Iteration $i -----------------
        echo ./build/heat_stencil_1D_seq $N
        ./build/heat_stencil_1D_seq $N
        echo
        echo mpiexec build/heat_stencil_1D_par $N
        mpiexec build/heat_stencil_1D_par $N
        echo
    done
done