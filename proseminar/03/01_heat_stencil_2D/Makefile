CC=gcc
MPICC=mpicc
CFLAGS=-O3 -std=gnu11 -Wall -Wextra -pedantic

.PHONY: all clean run process_data

all: heat_stencil_2D_seq heat_stencil_2D_par_blocking

heat_stencil_2D_seq: src/heat_stencil_2D_seq.c  
	@mkdir -p build
	@$(CC) $(CFLAGS) $< -o build/heat_stencil_2D_seq

heat_stencil_2D_par_blocking: src/heat_stencil_2D_par_blocking.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_2D_par_blocking -lm

clean:
	rm -rf build

clean_all:
	rm -rf build output

run:
	rm -rf build output
	./run.sh

process_data:
	module load Anaconda3/2023.10/python-3.11.5-anaconda-2023.10 && python3 src/process_data.py output/measurements.csv