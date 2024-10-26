CC=gcc
MPICC=mpicc
CFLAGS=-O3 -std=gnu11 -Wall -Wextra -pedantic -mtune=native

.PHONY: all clean run process_data

all: heat_stencil_2D_seq heat_stencil_2D_par_blocking heat_stencil_2D_par_blockingSlices heat_stencil_2D_seq_modified_access  heat_stencil_2D_par_blockingSlices_modified_access

heat_stencil_2D_seq: src/heat_stencil_2D_seq.c
	@mkdir -p build
	@$(CC) $(CFLAGS) $< -o build/heat_stencil_2D_seq

heat_stencil_2D_seq_modified_access: src/heat_stencil_2D_seq_modified_access.c
	@mkdir -p build
	@$(CC) $(CFLAGS) $< -o build/heat_stencil_2D_seq_modified_access

# heat_stencil_2D_seqFFT: src/heat_stencil_2D_seqFFT.c src/FFT_Conv2D.o
# 	@mkdir -p build
# 	@$(CC) $(CFLAGS) $< -o build/heat_stencil_2D_seqFFT

heat_stencil_2D_par_blocking: src/heat_stencil_2D_par_blocking.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_2D_par_blocking -lm

heat_stencil_2D_par_blockingSlices: src/heat_stencil_2D_par_blockingSlices.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_2D_par_blockingSlices

heat_stencil_2D_par_blockingSlices_modified_access: src/heat_stencil_2D_par_blockingSlices_modified_access.c
	@mkdir -p build
	@$(MPICC) $(CFLAGS) $< -o build/heat_stencil_2D_par_blockingSlices_modified_access

clean:
	rm -rf build

clean_all:
	rm -rf build output

run:
	./run.sh

process_data:
	module load Anaconda3/2023.10/python-3.11.5-anaconda-2023.10 && python3 src/process_data.py output/measurements.csv
