CC=gcc
MPICC=mpicc
CFLAGS=-O2 -std=gnu11 -Wall -Wextra -pedantic

.PHONY: all clean run plot

all: heat_stencil_1D_seq heat_stencil_1D_par_complex heat_stencil_1D_par_simple

heat_stencil_1D_seq: src/heat_stencil_1D_seq.c
	@mkdir -p build
	@$(CC) $(CFLAGS) $< -o build/heat_stencil_1D_seq

heat_stencil_1D_par_complex: src/heat_stencil_1D_par_complex.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_1D_par_complex

heat_stencil_1D_par_simple: src/heat_stencil_1D_par_simple.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_1D_par_simple

clean:
	rm -rf build output

run:
	rm -rf build output
	sbatch job.sh

plot:
	module load Anaconda3/2023.10/python-3.11.5-anaconda-2023.10 && python3 src/plot.py output/measurements.csv