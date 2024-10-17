#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <mpi.h>
#include <math.h>

#define RESOLUTION_WIDTH  192
#define RESOLUTION_HEIGHT 192

// ---------- VECTOR UTILITIES ----------
typedef double value_t;
typedef value_t *Vector;
Vector createVector(int N);
void releaseVector(Vector m);

// ---------- MATRIX UTILITIES ----------
int calc_index(int i, int j, int N);
int calc_index_supermatrix(int rank, int N, int ranks_per_row, int rows_per_rank, int cols_per_rank, int index_submatrix);
typedef value_t *Matrix;
Matrix createMatrix(int N, int M);
void calc_rank_factors(int numRanks, int *ranks_per_row, int *ranks_per_col);
void releaseMatrix(Matrix m);
void printTemperature(double *m, int N, int M);

// ---------- MEASUREMENT UTILITIES ----------
#define FOLDER "output"
#define FILENAME "measurements.csv"
#define TYPE "par_block"
void timings_to_csv(int problem_size, double time, int numRanks);

// ---------- SIMULATION CODE ----------
int main(int argc, char **argv) {
    // 'parsing' optional input parameter = problem size
    int N = 200;
    if (argc > 1) {
        N = atoi(argv[1]);
    }
    int T = N * 100;

    // ---------- START CLOCK ----------
    clock_t start = clock();

    // ---------- MPI INIT ----------
    MPI_Init(&argc, &argv);
    int myRank, numRanks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    if (N % numRanks) {
        printf("Configuration not possible: N=%d %% ranks=%d != 0\n", N, numRanks);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // ---------- SETUP ----------
    // calculate rank factors: ranks_per_row x ranks_per_col = numRanks
    int ranks_per_row = 0;
    int ranks_per_col = 0;
    calc_rank_factors(numRanks, &ranks_per_row, &ranks_per_col);
    printf("ranks per row %d and col %d\n", ranks_per_row, ranks_per_col);
    int rows_per_rank = N / ranks_per_col;
    int cols_per_rank = N / ranks_per_row;

    if (myRank == 0) {
        printf("Computing heat-distribution for room size %dX%d for T=%d timesteps\n", N, N, T);
        printf("Room is distributed across %d processes, each handling a sub-room of size %dX%d\n", numRanks, rows_per_rank, cols_per_rank);
    }

    // create a matrix for storing temperature fields
    Matrix A = NULL;
    // create a second matrix for the computation
    Matrix B = NULL;
    // create a third matrix for print/verification aggregation
    Matrix printMat = NULL;

    if (myRank == 0) {
        printMat = createMatrix(N, N);
        for (int i = 0; i < N*N; i++) {
           printMat[i] = 273;
        }
        // for(int i = 0; i<rows_per_rank*cols_per_rank; i++)

        printf("ranks per row %d \n ranks per col %d \n cols per rank %d \n rows per ranlk %d\n", ranks_per_row, ranks_per_col, cols_per_rank, rows_per_rank);

        for(int i=0; i<numRanks; i++)
        printMat[(1%ranks_per_row)*rows_per_rank+1/ranks_per_row*cols_per_rank*N] = 273+60; 
    }

    A = createMatrix(rows_per_rank, cols_per_rank);
    B = createMatrix(rows_per_rank, cols_per_rank);
    

    // set up initial conditions in A
    for (int i = 0; i < rows_per_rank*cols_per_rank; i++) {
        A[i] = 273; // temperature is 0° C everywhere (273 K)
        B[i] = 273; // this is not required TODO remove
    }

    // and there is a heat source
    int source_i = N / 4;
    int source_j = 190;

    int block_row = source_i / rows_per_rank;
    int block_col = source_j / cols_per_rank;
    int rank_with_source = block_row * ranks_per_row + block_col;
    printf("heat source located in rank %d\n", rank_with_source);

    if (myRank == 0) {
        printMat[calc_index(source_i, source_j, N)] = 273 + 60;
    }
    if (myRank == rank_with_source) {
        A[calc_index(source_i % rows_per_rank, source_j % cols_per_rank, cols_per_rank)] = 273 + 60;
    }

    if (myRank == 0) {
        printf("Initial:\n");
        printTemperature(printMat, N, N);
        printf("\n");
    }
    // create ghost vectors for the rank communication
    Vector left_temps = createVector(rows_per_rank);
    Vector right_temps = createVector(rows_per_rank);
    Vector up_temps = createVector(cols_per_rank);
    Vector down_temps = createVector(cols_per_rank);

    // ---------- COMPUTE ----------
    // for each time step ..
    for (int t = 0; t < T; t++) {
        // communication between ranks to get temperatures of adjacent cells
        // TODO

        // .. we propagate the temperature
        for (int i = 0; i < rows_per_rank; i++) {
            for (int j = 0; j < cols_per_rank; j++) {
                if(myRank == rank_with_source && (i == (source_i % rows_per_rank)) && (j == (source_j % cols_per_rank))){
                    B[calc_index(i, j, cols_per_rank)] = A[calc_index(i, j, cols_per_rank)];
                    continue;
                }

                // get temperature at current position
                // value_t current_temp = A[calc_index(i, j, rows_per_rank)];

                // // get temperatures of adjacent cells
                // value_t left_temp = (j != 0) ? A[calc_index(i, j-1, rows_per_rank)] : left_temps[i];
                // value_t right_temp = (j != N-1) ? A[calc_index(i, j+1, rows_per_rank)] : right_temps[i];
                // value_t up_temp = (i != 0) ? A[calc_index(i-1, j, rows_per_rank)] : up_temps[j];
                // value_t down_temp = (i != N-1) ? A[calc_index(i+1, j, rows_per_rank)] : down_temps[j];

                // B[calc_index(i, j, rows_per_rank)] = current_temp + 1/8.f * (left_temp + right_temp + down_temp + up_temp + (-4 * current_temp));
                // B[calc_index(i, j, rows_per_rank)] = (left_temp + right_temp + 4*current_temp + down_temp + up_temp)/8;
                // B[calc_index(i, j, rows_per_rank)] = 0;
            }
        }

        // swap matrices (just pointers, not content)
        Matrix H = A;
        A = B;
        B = H;

        // every 1000 steps show intermediate step
        if (!(t % 1000)) {
            // Gather all data from all ranks to rank 0
            int elements_per_rank = rows_per_rank*cols_per_rank;
            for(int i = 0; i<elements_per_rank; i++){
                A[i] = 278+myRank*5;
            }


            if (myRank != 0){
                MPI_Send(A, elements_per_rank, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
            if (myRank == 0) {
                for(int i = 0; i < elements_per_rank; i++){
                    int index_in_printMat = calc_index_supermatrix(myRank, N, ranks_per_row, rows_per_rank, cols_per_rank, i);
                    printMat[index_in_printMat] = A[i]; 
                    // printf("index in printmat: %d comming from rank %d\n", index_in_printMat ,i);
                    // printf("wrote into %d from %d\n", i%cols_per_rank+myRank*elements_per_rank+(i/cols_per_rank)*N, i);
                }
                Matrix tmp = createMatrix(rows_per_rank, cols_per_rank);
                for (int i = 1; i<numRanks; i++){
                    MPI_Recv(tmp, elements_per_rank, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, NULL);
                    for(int j = 0; j < elements_per_rank; j++){
                        int index_in_printMat = calc_index_supermatrix(i, N, ranks_per_row, rows_per_rank, cols_per_rank, j);
                        printMat[index_in_printMat] = tmp[j]; 
                    }
                    
                }
                releaseMatrix(tmp);
            }


            if (myRank == 0) {
                printf("Step t=%d\n", t);
                printTemperature(printMat, N, N);
                printf("\n");
            }
        }
    }

    // cleanup resources of computation
    releaseMatrix(B);
    releaseVector(left_temps);
    releaseVector(right_temps);
    releaseVector(up_temps);
    releaseVector(down_temps);

    // ---------- FINAL RESULT ----------
    int success = 1;
    // Gather all data from all ranks to rank 0
    // for (int i = 0; i < N; i++) {
    //     MPI_Gather(&A[calc_index(i, 0, rows_per_rank)], rows_per_rank, MPI_DOUBLE, &A[calc_index(i, 0, rows_per_rank)], rows_per_rank, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // }
    if (myRank == 0)  {
        printf("Final:\n");
        printTemperature(printMat, N, N);
        printf("\n");

        // ---------- STOP CLOCK ----------
        clock_t end = clock();
        double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
        timings_to_csv(N, total_time, numRanks);

        // ---------- VERIFICATION ----------
        // simple verification if nowhere the heat is more then the heat source
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double temp = printMat[calc_index(i, j, N)];
                if (273 <= temp && temp <= 273 + 60) {
                    continue;
                }
                success = 0;
                break;
            }
        }
        printf("Verification: %s\n", (success) ? "OK" : "FAILED");
        printf("Wall Clock Time = %f seconds\n", total_time);
        releaseMatrix(printMat); //only rank 0 possesses this mat
    }

    // ---------- CLEANUP ----------
    releaseMatrix(A);
    
    MPI_Finalize();
    // return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
    return 0;
}

// ---------- UTILITIES ----------
int calc_index(int i, int j, int N) {
    return ((i) * (N) + (j));
}

// calculates the index in the supermatrix based on the index in the submatrix of a specific rank
int calc_index_supermatrix(int rank, int N, int ranks_per_row, int rows_per_rank, int cols_per_rank, int index_submatrix){
    return ((rank / ranks_per_row) * (N * rows_per_rank) + (rank % ranks_per_row) * cols_per_rank) + (index_submatrix % cols_per_rank) + (index_submatrix/cols_per_rank) * cols_per_rank + (index_submatrix/cols_per_rank) * (N-cols_per_rank);
}

Vector createVector(int N) {
  // create data and index vector
  Vector v = malloc(sizeof(value_t) * N);
  for (int i = 0; i < N; i++) {
    v[i] = 10;
  }
  return v;
}

void releaseVector(Vector m) { free(m); }

Matrix createMatrix(int N, int M) {
  // create data and index matrix
  return malloc(sizeof(value_t) * N * M);
}

void calc_rank_factors(int numRanks, int *ranks_per_row, int *ranks_per_col) {
    // find ranks_per_row x ranks_per_col == numRanks
    // best result would be sqrt(numRanks) x sqrt(numRanks) == numRanks
    // otherwise find almost evenly squared matrices
    for (*ranks_per_row = (int)sqrt(numRanks); *ranks_per_row >= 1; (*ranks_per_row)--) {
        if (numRanks % *ranks_per_row == 0) {
            *ranks_per_col = numRanks / *ranks_per_row;
            break;
        }
    }
}

void releaseMatrix(Matrix m) { free(m); }

void printTemperature(double *m, int N, int M) {
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
                    max_t = (max_t < m[calc_index(x, y, N)]) ? m[calc_index(x, y, N)] : max_t;
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

void timings_to_csv(int problem_size, double time, int numRanks) {
    FILE* fpt;
    int set_header = 0;
    char full_filepath[1024];
    sprintf(full_filepath, "%s/%s", FOLDER, FILENAME);
    if(access(FOLDER, F_OK) != 0) mkdir(FOLDER, 0755);
	if(access(full_filepath, F_OK) != 0) set_header = 1;
	fpt = fopen(full_filepath, "a+");
	if(set_header) fprintf(fpt, "Impl/Ranks,Problem Size,Time\n");
	fprintf(fpt, "%s/%d,%u,%.9f\n", TYPE, numRanks, problem_size, time);
	fclose(fpt);
}