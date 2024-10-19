#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <mpi.h>

typedef double value_t;
#define RESOLUTION 120

// -- vector utilities --
typedef value_t *Vector;
Vector create_vector(int N);
void release_vector(Vector m);
void print_temperature(Vector m, int N);

// -- measurment utilities --
#define FOLDER "output"
#define FILENAME "measurements.csv"
#define TYPE "parallel"
void data_to_csv(unsigned problem_size, double time, int num_ranks);

// -- simulation code ---
int main(int argc, char **argv) {
  clock_t start = clock();

  if(argc < 2){
		printf("Usage: %s <number of iterations>\n", argv[0]);
		return EXIT_FAILURE;
	}

  int N = atoi(argv[1]);
  int T = N * 500;

  MPI_Init(&argc, &argv);
  int my_rank, num_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  if (N % num_ranks) {
    if (my_rank == 0)
      printf("Configuration not possible: N=%d, ranks=%d\n", N, num_ranks);
    MPI_Finalize();
    return EXIT_FAILURE;
  }
  if (my_rank == 0)
    printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields per rank
  int vector_size_per_rank = N / num_ranks;
  Vector A = NULL;
  Vector B = NULL;
  if (my_rank == 0) {
    A = create_vector(N);
    B = create_vector(N);
  } else {
    A = create_vector(vector_size_per_rank);
    B = create_vector(vector_size_per_rank);
  }

  // set up initial conditions in A
  for (int i = 0; i < vector_size_per_rank; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source somewhere
  int source_x = N / 4;
  int source_y = 273 + 60;

  int rank_with_source = source_x / vector_size_per_rank;
  if (my_rank == 0) {
    A[source_x] = source_y;
  }
  if (my_rank == rank_with_source) {
    A[source_x % vector_size_per_rank] = source_y;
  }

  if (my_rank == 0) {
    printf("Initial:\t");
    print_temperature(A, N);
    printf("\n");
  }
  MPI_Request requests[4];
  // ---------- compute ----------

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // sync with other ranks
    value_t t_left, t_right;
    if (my_rank != 0 && my_rank != num_ranks-1){ //consider edge cases separatelly
      MPI_Isend(&A[vector_size_per_rank-1], 1, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &requests[0]); // request is NULLptr as we're waiting in the recv step for all
      MPI_Isend(&A[0], 1, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &requests[1]);

      MPI_Irecv(&t_left, 1, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &requests[2]);
      MPI_Irecv(&t_right, 1, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &requests[3]);
      
      MPI_Wait(&requests[2], MPI_STATUS_IGNORE);
      MPI_Wait(&requests[3], MPI_STATUS_IGNORE);
    }
    else if(my_rank == 0){ // edge case 1 first rank
      MPI_Isend(&A[vector_size_per_rank-1], 1, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &requests[0]);
      MPI_Irecv(&t_right, 1, MPI_DOUBLE, my_rank+1, 0, MPI_COMM_WORLD, &requests[1]);
      t_left = A[0];
    }
    else{ // edge case 1 last rank
      MPI_Isend(&A[0], 1, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &requests[0]);
      MPI_Irecv(&t_left, 1, MPI_DOUBLE, my_rank-1, 0, MPI_COMM_WORLD, &requests[1]);
      t_right = A[vector_size_per_rank-1];
    }
    // pre-calc syncing done
    MPI_Wait(&requests[0], MPI_STATUS_IGNORE);
    MPI_Wait(&requests[1], MPI_STATUS_IGNORE);


    // .. we propagate the temperature
    for (long long i = 0; i < vector_size_per_rank; i++) {
      // center stays constant (the heat is still on)
      if (my_rank == rank_with_source && i == (source_x%vector_size_per_rank)) {
        B[i] = A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? A[i - 1] : t_left;
      value_t tr = (i != vector_size_per_rank - 1) ? A[i + 1] : t_right;

      // compute new temperature at current position
      B[i] = tc + 0.2 * (tl + tr + (-2 * tc));
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 10000)) {
      // now we have to gather all data in rank 0
      MPI_Gather(A, vector_size_per_rank, MPI_DOUBLE, A, vector_size_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if (my_rank == 0){
        printf("Step t=%d:\t", t);
        print_temperature(A, N);
        printf("\n");
      }
    }
  }

  release_vector(B);
  int success = 1;
  // last sync
  MPI_Gather(A, vector_size_per_rank, MPI_DOUBLE, A, vector_size_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (my_rank == 0){
    // measure time
    clock_t end = clock();
    double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    data_to_csv(N, total_time, num_ranks);

    // ---------- check ----------

    printf("Final:\t\t");
    print_temperature(A, N);
    printf("\n");

    
    for (long long i = 0; i < N; i++) {
      value_t temp = A[i];
      if (273 <= temp && temp <= 273 + 60)
        continue;
      success = 0;
      break;
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
    printf("Wall Clock Time = %f seconds\n", total_time);
  }

  // ---------- cleanup ----------

  release_vector(A);

  MPI_Finalize();
  // done
  return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Vector create_vector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void release_vector(Vector m) { free(m); }

void print_temperature(Vector m, int N) {
  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 60;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int W = RESOLUTION;

  // step size in each dimension
  int sW = N / W;

  // room
  // left wall
  printf("X");
  // actual room
  for (int i = 0; i < W; i++) {
    // get max temperature in this tile
    value_t max_t = 0;
    for (int x = sW * i; x < sW * i + sW; x++) {
      max_t = (max_t < m[x]) ? m[x] : max_t;
    }
    value_t temp = max_t;

    // pick the 'color'
    int c = ((temp - min) / (max - min)) * numColors;
    c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

    // print the average temperature
    printf("%c", colors[c]);
  }
  // right wall
  printf("X");
}

void data_to_csv(unsigned problem_size, double time, int num_ranks) {
	FILE* fpt;
	int set_header = 0;
  char full_filepath[1024];
  sprintf(full_filepath, "%s/%s", FOLDER, FILENAME);
  if(access(FOLDER, F_OK) != 0) mkdir(FOLDER, 0755);
	if(access(full_filepath, F_OK) != 0) set_header = 1;
	fpt = fopen(full_filepath, "a+");
	if(set_header) fprintf(fpt, "Impl/Ranks,Problem Size,Time\n");
	fprintf(fpt, "par_simple/%d,%u,%.9f\n", num_ranks, problem_size, time);
	fclose(fpt);
}