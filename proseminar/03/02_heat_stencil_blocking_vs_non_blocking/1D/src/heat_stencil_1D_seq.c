#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

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
#define TYPE "seq"
void data_to_csv(unsigned problem_size, double time, int num_ranks);

// -- simulation code ---
int main(int argc, char **argv) {
  clock_t start = clock();
  // 'parsing' optional input parameter = problem size
  int N = 2000;
  if (argc > 1) {
    N = atoi(argv[1]);
  }
  int T = N * 500;
  printf("Computing heat-distribution for room size N=%d for T=%d timesteps\n", N, T);

  // ---------- setup ----------

  // create a buffer for storing temperature fields
  Vector A = create_vector(N);

  // set up initial conditions in A
  for (int i = 0; i < N; i++) {
    A[i] = 273; // temperature is 0° C everywhere (273 K)
  }

  // and there is a heat source in one corner
  int source_x = N / 4;
  A[source_x] = 273 + 60;

  printf("Initial:\t");
  print_temperature(A, N);
  printf("\n");

  // ---------- compute ----------

  // create a second buffer for the computation
  Vector B = create_vector(N);

  // for each time step ..
  for (int t = 0; t < T; t++) {
    // .. we propagate the temperature
    for (long long i = 0; i < N; i++) {
      // center stays constant (the heat is still on)
      if (i == source_x) {
        B[i] = A[i];
        continue;
      }

      // get temperature at current position
      value_t tc = A[i];

      // get temperatures of adjacent cells
      value_t tl = (i != 0) ? A[i - 1] : tc;
      value_t tr = (i != N - 1) ? A[i + 1] : tc;

      // compute new temperature at current position
      B[i] = tc + 0.2 * (tl + tr + (-2 * tc));
    }

    // swap matrices (just pointers, not content)
    Vector H = A;
    A = B;
    B = H;

    // show intermediate step
    if (!(t % 10000)) {
      printf("Step t=%d:\t", t);
      print_temperature(A, N);
      printf("\n");
    }
  }

  release_vector(B);

  // measure time
  clock_t end = clock();
  double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
	data_to_csv(N, total_time, 1);

  // ---------- check ----------

  printf("Final:\t\t");
  print_temperature(A, N);
  printf("\n");

  int success = 1;
  for (long long i = 0; i < N; i++) {
    value_t temp = A[i];
    if (273 <= temp && temp <= 273 + 60)
      continue;
    success = 0;
    break;
  }

  printf("Verification: %s\n", (success) ? "OK" : "FAILED");
  printf("Wall Clock Time = %f seconds\n", total_time);

  // ---------- cleanup ----------

  release_vector(A);

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
	fprintf(fpt, "%s/%d,%u,%.9f\n", TYPE, num_ranks, problem_size, time);
	fclose(fpt);
}