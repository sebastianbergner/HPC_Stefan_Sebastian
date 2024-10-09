#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_NUM_THREADS 16
#define N 500000000

int THREADS;

void* mc_pi_thread(void* args) {
	unsigned int* seedptr = (unsigned int*)malloc(sizeof(unsigned int));
	int thread_id = (int)args;
	*seedptr = (unsigned int)(thread_id * pthread_self());
	int* in_count = (int*)malloc(sizeof(int));
	*in_count = 0;
	int iterations_max = N / THREADS;
	for(int i = 0; i < iterations_max; i++) {
		double x = rand_r(seedptr) / (double)RAND_MAX;
		double y = rand_r(seedptr) / (double)RAND_MAX;
		if(x * x + y * y <= 1) {
			(*in_count) += 1;
		}
	}

	if(thread_id == 0) {
		int rest = N % THREADS;
		for(int i = 0; i < rest; ++i) {
			double x = rand_r(seedptr) / (double)RAND_MAX;
			double y = rand_r(seedptr) / (double)RAND_MAX;
			if(x * x + y * y <= 1) {
				(*in_count)++;
			}
		}
	}
	free(seedptr);
	pthread_exit((void*)in_count);
}

int main(int argc, char* argv[]) {
	if(argc < 2 || atoi(argv[1]) > 16) {
		fprintf(stderr, "./monte_carlo_parallel <number of threads [MAX 16]>");
		return 1;
	}
	THREADS = atoi(argv[1]);
	pthread_t thread_arr[MAX_NUM_THREADS];

	for(int i = 0; i < THREADS; i++) {
		if(pthread_create(&(thread_arr[i]), NULL, mc_pi_thread, (void*)i)) {
			printf("Error creating thread %i", i);
			return -1;
		}
	}
	void* outval;
	int sum = 0;
	for(int i = 0; i < THREADS; i++) {
		pthread_join((thread_arr[i]), &outval);
		sum += *(int*)outval;
		free(outval);
	}
	printf("calculated approx pi with %i points and %i threads, result is: %f\n", N, THREADS,
	       4.0 * sum / N);
	return 0;
}