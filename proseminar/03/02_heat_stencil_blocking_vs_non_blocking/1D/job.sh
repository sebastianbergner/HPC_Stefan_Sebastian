#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps3.2_1D
#SBATCH --output=output/output.log
#SBATCH --time=1:00:00
# #SBATCH --mail-user=stefan.r.wagner@student.uibk.ac.at,sebastian.bergner@student.uibk.ac.at
# #SBATCH --mail-user=stefan.r.wagner@student.uibk.ac.at
# #SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --ntasks=96

#SBATCH --exclusive

# load all modules
module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v # if compiling has to be done here
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55 

make --no-print-directory

REP=5
PROBLEMS="768 1536 3072 6144"
RANKS="2 6 12 24 48 96"

for i in $(seq 1 $REP); do
    echo ----------------- Iteration $i -----------------
    for N in $PROBLEMS; do
        echo ./build/heat_stencil_1D_seq $N
        ./build/heat_stencil_1D_seq $N
        echo
        for P in $RANKS; do
            echo mpiexec -n $P build/heat_stencil_1D_par_blocking $N
            mpiexec -n $P build/heat_stencil_1D_par_blocking $N
            echo
            echo mpiexec -n $P build/heat_stencil_1D_par_non_blocking $N
            mpiexec -n $P build/heat_stencil_1D_par_non_blocking $N
            echo
        done
    done
done