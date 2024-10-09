#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>


void timigs_to_csv(char* filename, unsigned problemsize, double time) {
	FILE* fpt;
	int set_header = 0;
	if(access(filename, F_OK) != 0) set_header = 1;
	fpt = fopen(filename, "a+");
	if(set_header) fprintf(fpt, "Type,Problem size,Time\n");
	fprintf(fpt, "%s,%u,%.9f\n", "parallel", problemsize, time);
	fclose(fpt);
}

int mc_pi(int S) {
	int in_count = 0;
	for(unsigned i = 0; i < S; ++i) {
		const double x = rand() / (double)RAND_MAX;
		const double y = rand() / (double)RAND_MAX;
		if(x * x + y * y <= 1.f) {
			in_count++;
		}
	}
	return in_count; // only return the number of elements inside
}

int main(int argc, char* argv[]) {
	if(argc < 2){
		printf("usage: %s <number of iterations>\n", argv[0]);
		exit(-1);
	}
	int N = strtoul(argv[1], NULL, 10);
	clock_t start = clock();
	MPI_Init(&argc, &argv); //start mpi
	// let every process work on their problem (we won't instruct them from the root node 0 only gather the data)
	int myRank, numProcesses;
	int insideLocal = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
	srand(time(NULL)*myRank);

	int calculations_per_process = N/numProcesses;
	insideLocal = mc_pi(calculations_per_process);
	
	// here https://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/ was quite helpful
	int insideGlobal = 0;
	//MPI_Reduce( void* send_data, void* recv_data, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm communicator)
	MPI_Reduce(&insideLocal, &insideGlobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (myRank == 0){
		printf("inside local %d inside global %d\n\n", insideLocal, insideGlobal);
		double result = 4.f * insideGlobal / N;
		printf("pi approx with %i points: %f\n", N, result);
		clock_t end = clock();
		double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
		printf("took %f seconds\n", time_taken);
		timigs_to_csv("./measurements.csv", N, time_taken);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}