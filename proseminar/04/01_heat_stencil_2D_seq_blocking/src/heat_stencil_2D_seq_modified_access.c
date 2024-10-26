#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>

// ---------- PRINTING UTILITIES ----------
#define RESOLUTION_WIDTH  48
#define RESOLUTION_HEIGHT 48

// ---------- MATRIX UTILITIES ----------
#define IND(i, j) ((i) * (N) + (j))
typedef double value_t;
typedef value_t *Matrix;
Matrix create_matrix(int N, int M);
void release_matrix(Matrix m);
void print_temperature(double *m, int N, int M);

// ---------- MEASUREMENT UTILITIES ----------
#define FOLDER "output"
#define FILENAME "measurements.csv"
#define TYPE "seq_modaccess"
void data_to_csv(int problem_size, double time, int num_ranks);

// ---------- SIMULATION CODE ----------
int main(int argc, char **argv) {
    // 'parsing' optional input parameter = problem size
    int N = 200;
    if (argc > 1) {
        N = atoi(argv[1]);
    }
    int T = N * 100;

    

    // ---------- SETUP ----------
    // create a matrix for storing temperature fields
    Matrix A = create_matrix(N, N);
    // create a second matrix for the computation
    Matrix B = create_matrix(N, N);
    // set up initial conditions in A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[IND(i,j)] = 273; // temperature is 0° C everywhere (273 K)
        }
    }

    // and there is a heat source
    int source_x = N / 4;
    int source_y = N / 4;
    A[IND(source_x,source_y)] = 273 + 60;

    printf("Computing heat-distribution for room size %dX%d for T=%d timesteps\n", N, N, T);
    printf("Initial:\n");
    print_temperature(A, N, N);
    printf("\n");
    // ---------- START CLOCK ----------
    clock_t start = clock();
    // ---------- COMPUTE ----------
    // for each time step ..
    for (int t = 0; t < T; t++) {
        // .. we propagate the temperature
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if((i == source_x) && (j == source_y)){
                    B[IND(i, j)]=A[IND(i, j)];
                    continue;
                }

                // get temperature at current position
                int index_current = IND(i,j);
                value_t current_temp = A[index_current];

                // get temperatures of adjacent cells
                value_t left_temp = (j != 0) ? A[index_current-1] : current_temp;
                value_t right_temp = (j != N-1) ? A[index_current+1] : current_temp;
                value_t up_temp = (i != 0) ? A[index_current-N] : current_temp;
                value_t down_temp = (i != N-1) ? A[index_current+N] : current_temp;

                B[index_current] = current_temp + 1/8.f * (left_temp + right_temp + down_temp + up_temp + (-4 * current_temp));
            }
        }

        // swap matrices (just pointers, not content)
        Matrix H = A;
        A = B;
        B = H;

        // every 10000 steps show intermediate step
        if (!(t % 10000)) {
            printf("Step t=%d\n", t);
            print_temperature(A, N, N);
            printf("\n");
        }
    }

    

    // ---------- FINAL RESULT ----------
    printf("Final:\n");
    print_temperature(A, N, N);
    printf("\n");

    // ---------- STOP CLOCK ----------
    clock_t end = clock();
    double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    data_to_csv(N, total_time, 1);
    release_matrix(B);
    // ---------- VERIFICATION ----------
    // simple verification if nowhere the heat is more then the heat source
    int success = 1;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double temp = A[IND(i,j)];
            if (273 <= temp && temp <= 273 + 60) {
                continue;
            }
            success = 0;
            break;
        }
    }

    printf("Verification: %s\n", (success) ? "OK" : "FAILED");
    printf("Wall Clock Time = %f seconds\n", total_time);

    // ---------- CLEANUP ----------
    release_matrix(A);
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

// ---------- UTILITIES ----------
Matrix create_matrix(int N, int M) {
  // create data and index matrix
  return malloc(sizeof(value_t) * N * M);
}

void release_matrix(Matrix m) { free(m); }

void print_temperature(double *m, int N, int M) {
    const char *colors = " .-:=+*^X#%@";
    const int numColors = 12;

    // boundaries for temperature (for simplicity hard-coded)
    const double max = 273 + 60;
    const double min = 273 + 0;

    // set the 'render' resolution
    int W = RESOLUTION_WIDTH;
    int H = RESOLUTION_HEIGHT;

    // step size in each dimension
    int sW = N / W;
    int sH = M / H;

    // upper wall
    printf("\t");
    for (int u = 0; u < W + 2; u++) {
        printf("X");
    }
    printf("\n");
    // room
    for (int i = 0; i < H; i++) {
        // left wall
        printf("\tX");
        // actual room
        for (int j = 0; j < W; j++) {
            // get max temperature in this tile
            double max_t = 0;
            for (int x = sH * i; x < sH * i + sH; x++) {
                for (int y = sW * j; y < sW * j + sW; y++) {
                    max_t = (max_t < m[IND(x,y)]) ? m[IND(x,y)] : max_t;
                }
            }
            double temp = max_t;

            // pick the 'color'
            int c = ((temp - min) / (max - min)) * numColors;
            c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

            // print the average temperature
            printf("%c", colors[c]);
        }
        // right wall
        printf("X\n");
    }
    // lower wall
    printf("\t");
    for (int l = 0; l < W + 2; l++) {
        printf("X");
    }
    printf("\n");
}

void data_to_csv(int problem_size, double time, int num_ranks) {
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