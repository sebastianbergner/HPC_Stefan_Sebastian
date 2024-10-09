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
	fprintf(fpt, "%s,%u,%.9f\n", "serial", problemsize, time);
	fclose(fpt);
}

double mc_pi(unsigned int S) {
	int in_count = 0;
	for(unsigned i = 0; i < S; ++i) {
		const double x = rand() / (double)RAND_MAX;
		const double y = rand() / (double)RAND_MAX;
		if(x * x + y * y <= 1.f) {
			in_count++;
		}
	}
	return 4.f * in_count / S;
}

int main(int argc, char* argv[]) {
	if(argc < 2){
		printf("usage: %s <number of iterations>\n", argv[0]);
		exit(-1);
	}
	unsigned int N = strtoul(argv[1], NULL, 10);
	clock_t start = clock();
	printf("pi approx with %i points: %f\n", N, mc_pi(N));
	clock_t end = clock();
	double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
	printf("took %f seconds\n", time_taken);
	timigs_to_csv("./measurements.csv", N, time_taken);
	return 0;
}