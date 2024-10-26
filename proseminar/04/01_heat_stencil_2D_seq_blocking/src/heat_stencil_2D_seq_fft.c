#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>

// ---------- PRINTING UTILITIES ----------
#define RESOLUTION_WIDTH  20
#define RESOLUTION_HEIGHT 20

// ---------- MATRIX UTILITIES ----------
#define IND(i, j) ((i) * (N) + (j))
#define INDN(i, j, N) ((i) * (N) + (j))
typedef double value_t;
typedef value_t *Matrix;
Matrix create_matrix(int N, int M);
void release_matrix(Matrix m);
void print_temperature(double *m, int N, int M);

// ---------- MEASUREMENT UTILITIES ----------
#define FOLDER "output"
#define FILENAME "measurements.csv"
#define TYPE "seq_fft"
void data_to_csv(int problem_size, double time, int num_ranks);




// ------------- FFT UTILITIES -------------
typedef struct {
  float r;
  float i;
} complex;
static complex ctmp;

#define C_SWAP(a, b)                                                           \
  {                                                                            \
    ctmp = (a);                                                                \
    (a) = (b);                                                                 \
    (b) = ctmp;                                                                \
  }

void c_fft1d(complex *r, int n, int isign) {
  int m, i, i1, j, k, i2, l, l1, l2;
  float c1, c2, z;
  complex t, u;

  if (isign == 0)
    return;

  /* Do the bit reversal */
  i2 = n >> 1;
  j = 0;
  for (i = 0; i < n - 1; i++) {
    if (i < j)
      C_SWAP(r[i], r[j]);
    k = i2;
    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }

  /* m = (int) log2((double)n); */
  for (i = n, m = 0; i > 1; m++, i /= 2)
    ;

  /* Compute the FFT */
  c1 = -1.0;
  c2 = 0.0;
  l2 = 1;
  for (l = 0; l < m; l++) {
    l1 = l2;
    l2 <<= 1;
    u.r = 1.0;
    u.i = 0.0;
    for (j = 0; j < l1; j++) {
      for (i = j; i < n; i += l2) {
        i1 = i + l1;

        /* t = u * r[i1] */
        t.r = u.r * r[i1].r - u.i * r[i1].i;
        t.i = u.r * r[i1].i + u.i * r[i1].r;

        /* r[i1] = r[i] - t */
        r[i1].r = r[i].r - t.r;
        r[i1].i = r[i].i - t.i;

        /* r[i] = r[i] + t */
        r[i].r += t.r;
        r[i].i += t.i;
      }
      z = u.r * c1 - u.i * c2;

      u.i = u.r * c2 + u.i * c1;
      u.r = z;
    }
    c2 = sqrt((1.0 - c1) / 2.0);
    if (isign == -1) /* FWD FFT */
      c2 = -c2;
    c1 = sqrt((1.0 + c1) / 2.0);
  }

  /* Scaling for inverse transform */
  if (isign == 1) { /* IFFT*/
    for (i = 0; i < n; i++) {
      r[i].r /= n;
      r[i].i /= n;
    }
  }
}

void fft2d(complex* data, int isign, int N) {

  int i, j;
  complex *vec = (complex *)malloc(N * sizeof(complex));

  for (j = 0; j < N; j++) {
    for (i = 0; i < N; i++) {
      vec[i] = data[INDN(i, j, N)];
    }
    c_fft1d(vec, N, isign);
    for (i = 0; i < N; i++) {
      data[INDN(i, j, N)] = vec[i];
    }
  }

  free(vec);

  vec = (complex *)malloc(N * sizeof(complex));

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      vec[j] = data[INDN(i, j, N)];
    }
    c_fft1d(vec, N, isign);
    for (j = 0; j < N; j++) {
      data[INDN(i, j, N)] = vec[j];
    }
  }

  free(vec);
}

void mmpoint(complex* data1, complex* data2, complex* data3, int N) {
  int i, j;

  double real, imag = 0.f;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      data3[INDN(i, j, N)].r = (data1[INDN(i, j, N)].r * data2[INDN(i, j, N)].r) - (data1[INDN(i, j, N)].i * data2[INDN(i, j, N)].i);
      data3[INDN(i, j, N)].i = (data1[INDN(i, j, N)].r * data2[INDN(i, j, N)].i) + (data1[INDN(i, j, N)].i * data2[INDN(i, j, N)].r);
    }
  }
}

void double_to_complex(double* input, complex* output, const int n){
  for(int i = 0; i < n*n; i++){
    output[i].r = input[i];
    output[i].i = 0.0f;
  }
}

int next_power_of_2(int n) {
    int power = 1;
    while (power < n) {
        power *= 2;
    }
    return power;
}

void zero_pad(double* input, double* padded, int rows, int cols, int new_cols, int new_rows) {
    // Initialize all elements in padded to zero
    for (int i = 0; i < new_rows * new_cols; i++) {
        padded[i] = 0.0;
    }
    
    // Copy elements from input to padded within the bounds of rows and cols
    int pad_width_row = (new_rows-rows)/2;
    int pad_width_col = (new_cols-cols)/2;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            padded[INDN(i+pad_width_row, j+pad_width_col, new_rows)] = input[INDN(i, j, rows)];
        }
    }
}

void mirror_pad(double* input, double* padded, int rows, int cols, int new_cols, int new_rows) {
    int pad_width_row = (new_rows - rows) / 2;
    int pad_width_col = (new_cols - cols) / 2;

    // Initialize the padded array by mirroring the input content
    for (int i = 0; i < new_rows; i++) {
        for (int j = 0; j < new_cols; j++) {
            int src_i = i - pad_width_row;
            int src_j = j - pad_width_col;

            // Mirror the row index if out of bounds
            if (src_i < 0) src_i = -src_i;
            if (src_i >= rows) src_i = 2 * rows - src_i - 2;

            // Mirror the column index if out of bounds
            if (src_j < 0) src_j = -src_j;
            if (src_j >= cols) src_j = 2 * cols - src_j - 2;

            // Set the padded element
            padded[INDN(i, j, new_rows)] = input[INDN(src_i, src_j, rows)];
        }
    }
}


void printcomplexmat(complex* mat, int size){
for (int i = 0; i<size; i++){
    for (int j = 0; j<size; j++){
      printf("%lf ", mat[i * size+j].r);
    }
    printf("\n");
  }
}

void CircleShift(double* output, const double* input, int rows, int cols, int yshift, int xshift){
	for (int r = 0; r < rows; r++) {
		int newR = (r + yshift) % rows;
		if (newR < 0) newR = rows + newR;

		for (int c = 0; c < cols; c++) {
			int newC = (c + xshift) % cols;
			if (newC < 0) newC = cols + newC;

			output[newR * cols + newC] = input[r * cols + c];
		}
	}
}

void fft_shift(complex* data, int N) {
    int half_N = N / 2;
    complex temp;

    for (int i = 0; i < half_N; i++) {
        for (int j = 0; j < half_N; j++) {
            // Swap top-left with bottom-right
            temp = data[INDN(i, j, N)];
            data[INDN(i, j, N)] = data[INDN(i + half_N, j + half_N, N)];
            data[INDN(i + half_N, j + half_N, N)] = temp;

            // Swap top-right with bottom-left
            temp = data[INDN(i, j + half_N, N)];
            data[INDN(i, j + half_N, N)] = data[INDN(i + half_N, j, N)];
            data[INDN(i + half_N, j, N)] = temp;
        }
    }
}
// the results inputed here will be located in the left corner (remember zero padding) 
void complex_to_double(complex* input, double* output, int n_in, int n_out){
  for(int i = 0; i < n_out; i++){
    for(int j = 0; j< n_out; j++){
        output[INDN(i,j,n_out)] = input[INDN(i,j,n_in)].r;
    }
  }
}


void printmat(double* mat, int size, int stride){
for (int i = 0; i<size; i++){
    for (int j = 0; j<size; j++){
      printf("%lf ", mat[i * (size+stride)+j]);
    }
    printf("\n");
  }
}








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

    int optimal_fft_mat_size = next_power_of_2(N);
    // ---------- SETUP ----------
    // create a matrix for storing temperature fields
    Matrix A = create_matrix(N, N);
    // create a second matrix for the computation
    Matrix B = create_matrix(N, N);
    // set up initial conditions in A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[IND(i,j)] = 0; // temperature is 0° C everywhere (273 K)
            B[IND(i,j)] = 0;
        }
    }

    Matrix kernel_helper = create_matrix(3,3);
    for(int i = 0; i<3*3; i++){
        if (i == 1 || i == 3 || i == 5 || i == 7){
            kernel_helper[i] = 1/8.f;
        }
        else if(i == 4){
            kernel_helper[i] = .5f;
        }
        else{
            kernel_helper[i] = 0;
        }
    }

    Matrix A_helper = create_matrix(optimal_fft_mat_size,optimal_fft_mat_size); // for padding
    Matrix kernel = create_matrix(optimal_fft_mat_size,optimal_fft_mat_size); // not the most beautiful but hope it works
    zero_pad(kernel_helper, kernel, 3, 3, optimal_fft_mat_size,optimal_fft_mat_size); 
    
    complex* complex_kernel = (complex*)malloc(optimal_fft_mat_size*optimal_fft_mat_size*sizeof(complex));
    complex* complex_A = (complex*)malloc(optimal_fft_mat_size*optimal_fft_mat_size*sizeof(complex));
    complex* complex_B = (complex*)malloc(optimal_fft_mat_size*optimal_fft_mat_size*sizeof(complex));

    double_to_complex(kernel, complex_kernel, optimal_fft_mat_size);
    fft2d(complex_kernel, -1, optimal_fft_mat_size);

    // and there is a heat source
    int source_x = N / 4;
    int source_y = N / 4;
    A[IND(source_x,source_y)] = 273 + 60;

    printf("Computing heat-distribution for room size %dX%d for T=%d timesteps\n", N, N, T);
    printf("Initial:\n");
    print_temperature(A, N, N);
    printf("\n");

    // ---------- COMPUTE ----------
    // for each time step ..
    for (int t = 0; t < T; t++) {
        // .. we propagate the temperature
        // for (int i = 0; i < N; i++) {
        //     for (int j = 0; j < N; j++) {
        //         if((i == source_x) && (j == source_y)){
        //             B[IND(i, j)]=A[IND(i, j)];
        //             continue;
        //         }

        //         // get temperature at current position
        //         value_t current_temp = A[IND(i,j)];

        //         // get temperatures of adjacent cells
        //         value_t left_temp = (j != 0) ? A[IND(i, j-1)] : current_temp;
        //         value_t right_temp = (j != N-1) ? A[IND(i, j+1)] : current_temp;
        //         value_t up_temp = (i != 0) ? A[IND(i-1, j)] : current_temp;
        //         value_t down_temp = (i != N-1) ? A[IND(i+1, j)] : current_temp;

        //         B[IND(i,j)] = current_temp + 1/8.f * (left_temp + right_temp + down_temp + up_temp + (-4 * current_temp));
        //     }
        // }


        // fun starts here :) or not...

        zero_pad(A, A_helper, N, N, optimal_fft_mat_size, optimal_fft_mat_size);
        double_to_complex(A_helper, complex_A, optimal_fft_mat_size);
        
        fft2d(complex_A, -1, optimal_fft_mat_size);
        fft_shift(complex_A, optimal_fft_mat_size);
        mmpoint(complex_kernel, complex_A, complex_B, optimal_fft_mat_size);

        fft2d(complex_B, 1, optimal_fft_mat_size);
        
        fft_shift(complex_B, optimal_fft_mat_size);
        // printcomplexmat(complex_B, optimal_fft_mat_size);
        
        complex_to_double(complex_B, B, optimal_fft_mat_size, N);
        B[IND(source_x,source_y)] = 60;

        // swap matrices (just pointers, not content)
        Matrix H = A;
        A = B;
        B = H;

        // every 10000 steps show intermediate step
        if (!(t % 10000)) {
            for(int i = 0; i < N*N; i++)
                A[i]+=273;
            printf("Step t=%d\n", t);
            print_temperature(A, N, N);
            printf("\n");
            for(int i = 0; i < N*N; i++)
                A[i]-=273;
        }
    }

    release_matrix(B);

    // ---------- FINAL RESULT ----------
    printf("Final:\n");
    print_temperature(A, N, N);
    printf("\n");

    // ---------- STOP CLOCK ----------
    clock_t end = clock();
    double total_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    data_to_csv(N, total_time, 1);

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