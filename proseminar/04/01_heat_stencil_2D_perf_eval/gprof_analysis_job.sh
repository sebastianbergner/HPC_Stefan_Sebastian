#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps4.optim
#SBATCH --output=output/gprof_output.log
# SBATCH --open-mode=append
#SBATCH --time=00:10:00
#SBATCH --exclusive

#SBATCH --ntasks 6

module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

N=384

CFLAGS="-O1 -std=gnu11 -Wall -Wextra -pedantic"

mkdir -p build

gcc $CFLAGS -g src/heat_stencil_2D_seq.c -o build/heat_stencil_2D_seq_debug -pg
mpicc $CFLAGS src/heat_stencil_2D_par_blocking.c -o build/heat_stencil_2D_par_blocking -lm
mpicc $CFLAGS -g src/heat_stencil_2D_par_blocking.c -o build/heat_stencil_2D_par_blocking_debug -pg -lm
mpicc $CFLAGS src/heat_stencil_2D_par_non_blocking.c -o build/heat_stencil_2D_par_non_blocking -lm
mpicc $CFLAGS -g src/heat_stencil_2D_par_non_blocking.c -o build/heat_stencil_2D_par_non_blocking_debug -pg -lm

echo "compile done"

./build/heat_stencil_2D_seq_debug $N
mv gmon.out output/gmon_seq.out
mpiexec -n 1 build/heat_stencil_2D_par_blocking_debug $N : -n 5 build/heat_stencil_2D_par_blocking $N
mv gmon.out output/gmon_par_blocking6p.out
mpiexec -n 1 build/heat_stencil_2D_par_non_blocking_debug $N : -n 5 build/heat_stencil_2D_par_non_blocking $N
mv gmon.out output/gmon_par_non_blocking6p.out

gprof ./build/heat_stencil_2D_seq_debug output/gmon_seq.out > output/output_seq.txt
gprof ./build/heat_stencil_2D_par_blocking_debug output/gmon_par_blocking6p.out > output/output_par_blocking6p.txt
gprof ./build/heat_stencil_2D_par_non_blocking_debug output/gmon_par_non_blocking6p.out > output/output_par_non_blocking6p.txt