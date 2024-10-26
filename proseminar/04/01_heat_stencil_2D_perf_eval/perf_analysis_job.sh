#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps4.optim
#SBATCH --output=output/perf_output.log
# SBATCH --open-mode=append
#SBATCH --time=00:10:00
#SBATCH --exclusive

module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

N=384

CFLAGS="-O0 -std=gnu11 -Wall -Wextra -pedantic"

mkdir -p build

gcc $CFLAGS -g src/heat_stencil_2D_seq.c -o build/heat_stencil_2D_seq_debug -pg
mpicc $CFLAGS src/heat_stencil_2D_par_blocking.c -o build/heat_stencil_2D_par_blocking -lm
mpicc $CFLAGS -g src/heat_stencil_2D_par_blocking.c -o build/heat_stencil_2D_par_blocking_debug -pg -lm
mpicc $CFLAGS src/heat_stencil_2D_par_non_blocking.c -o build/heat_stencil_2D_par_non_blocking -lm
mpicc $CFLAGS -g src/heat_stencil_2D_par_non_blocking.c -o build/heat_stencil_2D_par_non_blocking_debug -pg -lm

echo "compile done"

perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads ./build/heat_stencil_2D_seq_debug $N
mpiexec -n 1 perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads build/heat_stencil_2D_par_blocking_debug $N : -n 5 build/heat_stencil_2D_par_blocking $N
mpiexec -n 1 perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads build/heat_stencil_2D_par_non_blocking_debug $N : -n 5 build/heat_stencil_2D_par_non_blocking $N

echo "perf done"