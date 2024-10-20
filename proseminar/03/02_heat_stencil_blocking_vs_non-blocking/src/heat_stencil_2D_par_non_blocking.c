#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <mpi.h>
#include <math.h>

// ---------- PRINTING UTILITIES ----------
#define RESOLUTION_WIDTH  48
#define RESOLUTION_HEIGHT 48

// ---------- MPI UTILITIES ----------
#define MPI_TAG 0

// ---------- VECTOR UTILITIES ----------
typedef double value_t;
typedef value_t *Vector;
Vector create_vector(int N);
void release_vector(Vector m);

// ---------- MATRIX UTILITIES ----------
int calc_index(int i, int j, int N);
int calc_index_supermatrix(int rank, int N, int ranks_per_row, int rows_per_rank, int cols_per_rank, int index_submatrix);
typedef value_t *Matrix;
Matrix create_matrix(int N, int M);
void calc_rank_factors(int num_ranks, int *ranks_per_row, int *ranks_per_col);
void release_matrix(Matrix m);
void gather_data_from_submatrices(Matrix submatrix, Matrix supermatrix, int my_rank, int num_ranks, int rows_per_rank, int cols_per_rank, int ranks_per_row, int N);
void print_temperature(double *m, int N, int M);

// ---------- MEASUREMENT UTILITIES ----------
#define FOLDER "output"
#define FILENAME "measurements.csv"
#define TYPE "par_non_block"
void data_to_csv(int problem_size, double time, int num_ranks);

// ---------- COMPUTATION UTILITIES ----------
void propagate_heat(int i_start, int i_max, int i_add, int j_start, int j_max, int j_add, Matrix A, Matrix B, int my_rank, int rank_with_source, int rows_per_rank, int cols_per_rank, int source_i, int source_j, int left, int right, int up, int down, Vector recv_left_temps, Vector recv_right_temps, Vector recv_up_temps, Vector recv_down_temps);

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
    int my_rank, num_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    if (N % num_ranks) {
        printf("Configuration not possible: N=%d %% ranks=%d != 0\n", N, num_ranks);
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // ---------- SETUP ----------
    // calculate rank factors: ranks_per_row x ranks_per_col = num_ranks
    int ranks_per_row = 0;
    int ranks_per_col = 0;
    calc_rank_factors(num_ranks, &ranks_per_row, &ranks_per_col);
    int rows_per_rank = N / ranks_per_col;
    int cols_per_rank = N / ranks_per_row;

    int submatrix_col_id = my_rank % ranks_per_row; 
    // int submatrix_row_id = my_rank / ranks_per_row; // not needed in non-blocking version

    if (my_rank == 0) {
        printf("Computing heat-distribution for room size %dX%d for T=%d timesteps\n", N, N, T);
        printf("Room is distributed across %d processes, each handling a sub-room of size %dX%d\n", num_ranks, rows_per_rank, cols_per_rank);
    }

    // create a matrix for storing temperature fields
    Matrix A = create_matrix(rows_per_rank, cols_per_rank);
    // create a second matrix for the computation
    Matrix B = create_matrix(rows_per_rank, cols_per_rank);
    // create a third matrix for print/verification aggregation
    Matrix print_matrix = NULL;
    if (my_rank == 0) {
        print_matrix = create_matrix(N, N);
        for (int i = 0; i < N*N; i++) {
           print_matrix[i] = 273;
        }
    }

    // set up initial conditions in A
    for (int i = 0; i < rows_per_rank*cols_per_rank; i++) {
        A[i] = 273; // temperature is 0° C everywhere (273 K)
    }

    // and there is a heat source
    int source_i = N / 4;
    int source_j = N / 4;

    int block_row = source_i / rows_per_rank;
    int block_col = source_j / cols_per_rank;
    int rank_with_source = block_row * ranks_per_row + block_col;

    if (my_rank == 0) {
        print_matrix[calc_index(source_i, source_j, N)] = 273 + 60;
    }
    if (my_rank == rank_with_source) {
        A[calc_index(source_i % rows_per_rank, source_j % cols_per_rank, cols_per_rank)] = 273 + 60;
    }

    if (my_rank == 0) {
        printf("Initial:\n");
        print_temperature(print_matrix, N, N);
        printf("\n");
    }
    
    // create send ghost vectors for the rank communication
    Vector send_up_temps = create_vector(cols_per_rank);
    Vector send_down_temps = create_vector(cols_per_rank);
    Vector send_left_temps = create_vector(rows_per_rank);
    Vector send_right_temps = create_vector(rows_per_rank);
    
    // create recv ghost vectors for the rank communication
    Vector recv_down_temps = create_vector(cols_per_rank);
    Vector recv_up_temps = create_vector(cols_per_rank);
    Vector recv_right_temps = create_vector(rows_per_rank);
    Vector recv_left_temps = create_vector(rows_per_rank);

    MPI_Request requests[8];

    // ---------- COMPUTE ----------
    // for each time step ..
    for (int t = 0; t < T; t++) {
        // gather send ghost vector data for all sides (4 loops because it reduces cache misses)
        for (int i = 0; i < cols_per_rank; i++){
            send_up_temps[i] = A[calc_index(0, i, cols_per_rank)];
        }
        for (int i = 0; i < cols_per_rank; i++){
            send_down_temps[i] = A[calc_index(rows_per_rank-1, i, cols_per_rank)];
        }
        for (int i = 0; i < rows_per_rank; i++){
            send_left_temps[i] = A[calc_index(i, 0, cols_per_rank)];
        }
        for (int i = 0; i < rows_per_rank; i++){
            send_right_temps[i] = A[calc_index(i, cols_per_rank-1, cols_per_rank)];
        }
        
        // determine neighbors
        int up    = (my_rank-ranks_per_row >= 0) ? my_rank-ranks_per_row : MPI_PROC_NULL;
        int down  = my_rank+ranks_per_row < num_ranks ? my_rank+ranks_per_row : MPI_PROC_NULL;
        int left  = submatrix_col_id ? my_rank-1 : MPI_PROC_NULL;
        int right = submatrix_col_id+1 < ranks_per_row ? my_rank+1 : MPI_PROC_NULL;

        // send ghost cell data to neighbors
        MPI_Isend(send_up_temps,    cols_per_rank, MPI_DOUBLE, up,    MPI_TAG, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(send_down_temps,  cols_per_rank, MPI_DOUBLE, down,  MPI_TAG, MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(send_left_temps,  rows_per_rank, MPI_DOUBLE, left,  MPI_TAG, MPI_COMM_WORLD, &requests[2]);
        MPI_Isend(send_right_temps, rows_per_rank, MPI_DOUBLE, right, MPI_TAG, MPI_COMM_WORLD, &requests[3]);
        
        // receive ghost cell data from neighbors
        MPI_Irecv(recv_up_temps,    cols_per_rank, MPI_DOUBLE, up,    MPI_TAG, MPI_COMM_WORLD, &requests[4]);
        MPI_Irecv(recv_down_temps,  cols_per_rank, MPI_DOUBLE, down,  MPI_TAG, MPI_COMM_WORLD, &requests[5]);
        MPI_Irecv(recv_left_temps,  rows_per_rank, MPI_DOUBLE, left,  MPI_TAG, MPI_COMM_WORLD, &requests[6]);
        MPI_Irecv(recv_right_temps, rows_per_rank, MPI_DOUBLE, right, MPI_TAG, MPI_COMM_WORLD, &requests[7]);

        // .. we propagate the temperature; do all possible calculations here
        int i_start = 1, i_max = (rows_per_rank-1), i_add = 1, j_start = 1, j_max = (cols_per_rank-1), j_add = 1;
        propagate_heat(i_start, i_max, i_add, j_start, j_max, j_add, A, B, my_rank, rank_with_source, rows_per_rank, cols_per_rank, source_i, source_j, left, right, up, down, recv_left_temps, recv_right_temps, recv_up_temps, recv_down_temps);

        // wait for all non-blocking communication to finish
        MPI_Waitall(8, requests, MPI_STATUS_IGNORE);

        // ... now propagate the heat for the up and down egde cases
        i_start = 0, i_max = rows_per_rank, i_add = (rows_per_rank-1), j_start = 0, j_max = cols_per_rank, j_add = 1;
        propagate_heat(i_start, i_max, i_add, j_start, j_max, j_add, A, B, my_rank, rank_with_source, rows_per_rank, cols_per_rank, source_i, source_j, left, right, up, down, recv_left_temps, recv_right_temps, recv_up_temps, recv_down_temps);

        // ... now propagate the heat for the left and right egde cases
        i_start = 0, i_max = rows_per_rank, i_add = 1, j_start = 0, j_max = cols_per_rank, j_add = (cols_per_rank-1);
        propagate_heat(i_start, i_max, i_add, j_start, j_max, j_add, A, B, my_rank, rank_with_source, rows_per_rank, cols_per_rank, source_i, source_j, left, right, up, down, recv_left_temps, recv_right_temps, recv_up_temps, recv_down_temps);

        // swap matrices (just pointers, not content)
        Matrix H = A;
        A = B;
        B = H;

        // every 10000 steps show intermediate step
        if (!(t % 10000)) {
            // Gather all data from all ranks to rank 0
            gather_data_from_submatrices(A, print_matrix, my_rank, num_ranks, rows_per_rank, cols_per_rank, ranks_per_row, N);

            if (my_rank == 0) {
                printf("Step t=%d\n", t);
                print_temperature(print_matrix, N, N);
                printf("\n");
            }
        }
    }

    // cleanup resources of computation
    release_matrix(B);
    release_vector(send_left_temps);
    release_vector(send_right_temps);
    release_vector(send_up_temps);
    release_vector(send_down_temps);

    release_vector(recv_right_temps);
    release_vector(recv_left_temps);
    release_vector(recv_down_temps);
    release_vector(recv_up_temps);

    // ---------- FINAL RESULT ----------
    int success = 1;

    // Gather all data from all ranks to rank 0
    gather_data_from_submatrices(A, print_matrix, my_rank, num_ranks, rows_per_rank, cols_per_rank, ranks_per_row, N);

    if (my_rank == 0)  {
        printf("Final:\n");
        print_temperature(print_matrix, N, N);
        printf("\n");

        // ---------- STOP CLOCK ----------
        clock_t end = clock();
        double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
        data_to_csv(N, total_time, num_ranks);

        // ---------- VERIFICATION ----------
        // simple verification if nowhere the heat is more then the heat source
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double temp = print_matrix[calc_index(i, j, N)];
                if (273 <= temp && temp <= 273 + 60) {
                    continue;
                }
                success = 0;
                break;
            }
        }
        printf("Verification: %s\n", (success) ? "OK" : "FAILED");
        printf("Wall Clock Time = %f seconds\n", total_time);
        release_matrix(print_matrix); //only rank 0 possesses this mat
    }

    // ---------- CLEANUP ----------
    release_matrix(A);
    
    MPI_Finalize();
    return (success) ? EXIT_SUCCESS : EXIT_FAILURE;
}

// ---------- UTILITIES ----------
int calc_index(int i, int j, int N) {
    return ((i) * (N) + (j));
}

// calculates the index in the supermatrix based on the index in the submatrix of a specific rank
int calc_index_supermatrix(int rank, int N, int ranks_per_row, int rows_per_rank, int cols_per_rank, int index_submatrix){
    return ((rank / ranks_per_row) * (N * rows_per_rank) + (rank % ranks_per_row) * cols_per_rank) + (index_submatrix % cols_per_rank) + (index_submatrix/cols_per_rank) * cols_per_rank + (index_submatrix/cols_per_rank) * (N-cols_per_rank);
}

Vector create_vector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void release_vector(Vector m) { free(m); }

Matrix create_matrix(int N, int M) {
  // create data and index matrix
  return malloc(sizeof(value_t) * N * M);
}

void calc_rank_factors(int num_ranks, int *ranks_per_row, int *ranks_per_col) {
    // find ranks_per_row x ranks_per_col == num_ranks
    // best result would be sqrt(num_ranks) x sqrt(num_ranks) == num_ranks
    // otherwise find almost evenly squared matrices
    for (*ranks_per_row = (int)sqrt(num_ranks); *ranks_per_row >= 1; (*ranks_per_row)--) {
        if (num_ranks % *ranks_per_row == 0) {
            *ranks_per_col = num_ranks / *ranks_per_row;
            break;
        }
    }
}

void release_matrix(Matrix m) { free(m); }

void gather_data_from_submatrices(Matrix submatrix, Matrix supermatrix, int my_rank, int num_ranks, int rows_per_rank, int cols_per_rank, int ranks_per_row, int N) {
    int elements_per_rank = rows_per_rank*cols_per_rank;

    if (my_rank != 0){
        MPI_Send(submatrix, elements_per_rank, MPI_DOUBLE, 0, MPI_TAG, MPI_COMM_WORLD);
    }
    if (my_rank == 0) {
        for(int i = 0; i < elements_per_rank; i++){
            int index_in_target = calc_index_supermatrix(my_rank, N, ranks_per_row, rows_per_rank, cols_per_rank, i);
            supermatrix[index_in_target] = submatrix[i]; 
        }
        Matrix tmp = create_matrix(rows_per_rank, cols_per_rank);
        for (int i = 1; i<num_ranks; i++){
            MPI_Recv(tmp, elements_per_rank, MPI_DOUBLE, i, MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j = 0; j < elements_per_rank; j++){
                int index_in_target = calc_index_supermatrix(i, N, ranks_per_row, rows_per_rank, cols_per_rank, j);
                supermatrix[index_in_target] = tmp[j]; 
            }
        }
        release_matrix(tmp);
    }
}

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

void propagate_heat(int i_start, int i_max, int i_add, int j_start, int j_max, int j_add, Matrix A, Matrix B, int my_rank, int rank_with_source, int rows_per_rank, int cols_per_rank, int source_i, int source_j, int left, int right, int up, int down, Vector recv_left_temps, Vector recv_right_temps, Vector recv_up_temps, Vector recv_down_temps) {
    for (int i = i_start; i < i_max; i += i_add) {
        for (int j = j_start; j < j_max; j += j_add) {
            if(my_rank == rank_with_source && (i == (source_i % rows_per_rank)) && (j == (source_j % cols_per_rank))){
                B[calc_index(i, j, cols_per_rank)] = A[calc_index(i, j, cols_per_rank)];
                continue;
            }

            // get temperature at current position
            value_t current_temp = A[calc_index(i, j, cols_per_rank)];

            // // get temperatures of adjacent cells
            value_t left_temp = (j != 0) ? A[calc_index(i, j-1, cols_per_rank)] : (left != MPI_PROC_NULL ? recv_left_temps[i] : current_temp);
            value_t right_temp = (j != cols_per_rank-1) ? A[calc_index(i, j+1, cols_per_rank)] : (right != MPI_PROC_NULL ? recv_right_temps[i] : current_temp);
            value_t up_temp = (i != 0) ? A[calc_index(i-1, j, cols_per_rank)] : (up != MPI_PROC_NULL ? recv_up_temps[j] : current_temp);
            value_t down_temp = (i != rows_per_rank-1) ? A[calc_index(i+1, j, cols_per_rank)] : (down != MPI_PROC_NULL ? recv_down_temps[j] : current_temp);
            
            B[calc_index(i, j, cols_per_rank)] = current_temp + 1/8.f * (left_temp + right_temp + down_temp + up_temp + (-4 * current_temp));
        }
    }
}