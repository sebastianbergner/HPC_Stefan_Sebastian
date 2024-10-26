#!/bin/bash
#SBATCH --partition=lva
#SBATCH --job-name ps4.optim
#SBATCH --output=output/output.log
# SBATCH --open-mode=append
#SBATCH --time=01:00:00
#SBATCH --mail-user=sebastian.bergner@student.uibk.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL
# stefan.r.wagner@student.uibk.ac.at,
#SBATCH --exclusive

module purge
module load gcc/12.2.0-gcc-8.5.0-p4pe45v
module load openmpi/3.1.6-gcc-12.2.0-d2gmn55

N=384

gcc -O0 -g src/heat_stencil_2D_seq.c -o build/heat_stencil_2D_seq -pg
gcc -O0 -g src/heat_stencil_2D_seq_modified_access.c -o build/heat_stencil_2D_seq_modified_access -pg
mpicc -O0 -g src/heat_stencil_2D_par_blocking.c -o build/heat_stencil_2D_par_blocking -pg -lm
mpicc -O0 -g src/heat_stencil_2D_par_blockingSlices.c -o build/heat_stencil_2D_par_blockingSlices -pg -lm
mpicc -O0 -g src/heat_stencil_2D_par_blockingSlices_modified_access.c -o build/heat_stencil_2D_par_blockingSlices_modified_access -pg -lm

echo "compile done"

./build/heat_stencil_2D_seq $N
mv gmon.out analysis/gmon_seq.out
./build/heat_stencil_2D_seq_modified_access $N
mv gmon.out analysis/gmon_seq_modified_access.out
mpiexec -n 6 build/heat_stencil_2D_par_blocking $N
mv gmon.out analysis/gmon_par_blocking6p.out
mpiexec -n 6 build/heat_stencil_2D_par_blockingSlices $N
mv gmon.out analysis/gmon_par_blocking_slices6p.out
mpiexec -n 6 build/heat_stencil_2D_par_blockingSlices_modified_access $N
mv gmon.out analysis/gmon_par_blockingSlices_modified_access6p.out

echo "gprof done"

perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads ./build/heat_stencil_2D_seq $N
perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads ./build/heat_stencil_2D_seq_modified_access $N
perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads mpiexec -n 6 build/heat_stencil_2D_par_blocking $N
perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads mpiexec -n 6 build/heat_stencil_2D_par_blockingSlices $N
perf stat -e L1-dcache-load-misses,L1-dcache-loads,LLC-load-misses,LLC-loads mpiexec -n 6 build/heat_stencil_2D_par_blockingSlices_modified_access $N

echo "perf done"

gcc -O3 src/heat_stencil_2D_seq.c -o build/heat_stencil_2D_seq
gcc -O3 src/heat_stencil_2D_seq_modified_access.c -o build/heat_stencil_2D_seq_modified_access
mpicc -O3 src/heat_stencil_2D_par_blocking.c -o build/heat_stencil_2D_par_blocking -lm
mpicc -O3 src/heat_stencil_2D_par_blockingSlices.c -o build/heat_stencil_2D_par_blockingSlices -lm
mpicc -O3 src/heat_stencil_2D_par_blockingSlices_modified_access.c -o build/heat_stencil_2D_par_blockingSlices_modified_access -lm

echo "compile 2 done"

./build/heat_stencil_2D_seq $N
./build/heat_stencil_2D_seq $N
mpiexec -n 6 build/heat_stencil_2D_par_blocking $N
mpiexec -n 6 build/heat_stencil_2D_par_blockingSlices $N
mpiexec -n 6 build/heat_stencil_2D_par_blockingSlices_modified_access $N

echo "execution 2 done"

mv ./output/measurements.csv ./analysis/measurements_analysis.csv
mv ./output/output.log ./analysis/output_analysis.log
