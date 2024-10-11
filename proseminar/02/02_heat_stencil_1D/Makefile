CC=gcc
MPICC=mpicc
CFLAGS=-O2 -std=gnu11 -Wall -Wextra -pedantic

.PHONY: all clean run plot

all: heat_stencil_1D_seq heat_stencil_1D_par

heat_stencil_1D_seq: src/heat_stencil_1D_seq.c
	@mkdir -p build
	@$(CC) $(CFLAGS) $< -o build/heat_stencil_1D_seq

heat_stencil_1D_par: src/heat_stencil_1D_par.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_1D_par

clean:
	rm -rf build output

run:
	rm -rf build output
	sbatch job.sh

plot:
	module load Anaconda3/2023.10/python-3.11.5-anaconda-2023.10 && python3 src/plot.py output/measurements.csv