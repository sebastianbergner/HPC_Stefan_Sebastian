#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>

#define INDN(i, j, N) ((i) * (N) + (j))

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
#define TYPE "seq_fft"
void data_to_csv(int problem_size, double time, int num_ranks);


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

// void c_fft1d(complex *r, int n, int isign) {
//   int m, i, i1, j, k, i2, l, l1, l2;
//   float c1, c2, z;
//   complex t, u;

//   if (isign == 0)
//     return;

//   /* Do the bit reversal */
//   i2 = n >> 1;
//   j = 0;
//   for (i = 0; i < n - 1; i++) {
//     if (i < j)
//       C_SWAP(r[i], r[j]);
//     k = i2;
//     while (k <= j) {
//       j -= k;
//       k >>= 1;
//     }
//     j += k;
//   }

//   /* m = (int) log2((double)n); */
//   for (i = n, m = 0; i > 1; m++, i /= 2)
//     ;

//   /* Compute the FFT */
//   c1 = -1.0;
//   c2 = 0.0;
//   l2 = 1;
//   for (l = 0; l < m; l++) {
//     l1 = l2;
//     l2 <<= 1;
//     u.r = 1.0;
//     u.i = 0.0;
//     for (j = 0; j < l1; j++) {
//       for (i = j; i < n; i += l2) {
//         i1 = i + l1;

//         /* t = u * r[i1] */
//         t.r = u.r * r[i1].r - u.i * r[i1].i;
//         t.i = u.r * r[i1].i + u.i * r[i1].r;

//         /* r[i1] = r[i] - t */
//         r[i1].r = r[i].r - t.r;
//         r[i1].i = r[i].i - t.i;

//         /* r[i] = r[i] + t */
//         r[i].r += t.r;
//         r[i].i += t.i;
//       }
//       z = u.r * c1 - u.i * c2;

//       u.i = u.r * c2 + u.i * c1;
//       u.r = z;
//     }
//     c2 = sqrt((1.0 - c1) / 2.0);
//     if (isign == -1) /* FWD FFT */
//       c2 = -c2;
//     c1 = sqrt((1.0 + c1) / 2.0);
//   }

//   /* Scaling for inverse transform */
//   if (isign == 1) { /* IFFT*/
//     for (i = 0; i < n; i++) {
//       r[i].r /= n;
//       r[i].i /= n;
//     }
//   }
// }

// void fft2d(complex* data, int isign, int N) {

//   int i, j;
//   complex *vec = (complex *)malloc(N * sizeof(complex));

//   for (j = 0; j < N; j++) {
//     for (i = 0; i < N; i++) {
//       vec[i] = data[INDN(i, j, N)];
//     }
//     c_fft1d(vec, N, isign);
//     for (i = 0; i < N; i++) {
//       data[INDN(i, j, N)] = vec[i];
//     }
//   }

//   free(vec);

//   vec = (complex *)malloc(N * sizeof(complex));

//   for (i = 0; i < N; i++) {
//     for (j = 0; j < N; j++) {
//       vec[j] = data[INDN(i, j, N)];
//     }
//     c_fft1d(vec, N, isign);
//     for (j = 0; j < N; j++) {
//       data[INDN(i, j, N)] = vec[j];
//     }
//   }

//   free(vec);
// }

#define M_PI 3.14159265358979323846
void fft(complex buf[], complex out[], int n, int step) {
    if (step < n) {
        fft(out, buf, n, step * 2);
        fft(out + step, buf + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step) {
            complex t = {cos(-M_PI * i / n), sin(-M_PI * i / n)};
            t.r *= out[i + step].r;
            t.i *= out[i + step].i;

            buf[i / 2].r = out[i].r + t.r;
            buf[i / 2].i = out[i].i + t.i;

            buf[(i + n) / 2].r = out[i].r - t.r;
            buf[(i + n) / 2].i = out[i].i - t.i;
        }
    }
}

void fft2d(complex matrix[], int rows, int cols) {
    complex *row = (complex *)malloc(cols * sizeof(complex));
    complex *col = (complex *)malloc(rows * sizeof(complex));

    if (row == NULL || col == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            row[c] = matrix[r * cols + c];
        }
        fft(row, row, cols, 1);
        for (int c = 0; c < cols; c++) {
            matrix[r * cols + c] = row[c];
        }
    }

    for (int c = 0; c < cols; c++) {
        for (int r = 0; r < rows; r++) {
            col[r] = matrix[r * cols + c];
        }
        fft(col, col, rows, 1);
        for (int r = 0; r < rows; r++) {
            matrix[r * cols + c] = col[r];
        }
    }

    free(row);
    free(col);
}

void ifft(complex buf[], complex out[], int n, int step) {
    if (step < n) {
        ifft(out, buf, n, step * 2);
        ifft(out + step, buf + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step) {
            complex t = {cos(M_PI * i / n), sin(M_PI * i / n)};
            t.r *= out[i + step].r;
            t.i *= out[i + step].i;

            buf[i / 2].r = out[i].r + t.r;
            buf[i / 2].i = out[i].i + t.i;

            buf[(i + n) / 2].r = out[i].r - t.r;
            buf[(i + n) / 2].i = out[i].i - t.i;
        }
    }
}

void ifft2d(complex matrix[], int rows, int cols) {
    complex *row = (complex *)malloc(cols * sizeof(complex));
    complex *col = (complex *)malloc(rows * sizeof(complex));

    if (row == NULL || col == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            row[c] = matrix[r * cols + c];
        }
        ifft(row, row, cols, 1);
        for (int c = 0; c < cols; c++) {
            matrix[r * cols + c] = row[c];
        }
    }

    for (int c = 0; c < cols; c++) {
        for (int r = 0; r < rows; r++) {
            col[r] = matrix[r * cols + c];
        }
        ifft(col, col, rows, 1);
        for (int r = 0; r < rows; r++) {
            matrix[r * cols + c] = col[r];
        }
    }

    free(row);
    free(col);
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

void to_complex(double* input, complex* output, const int n){
  for(int i = 0; i < n*n; i++){
    // printf("writing %")
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

void remove_zero_padding(double* padded, double* output, int rows, int cols, int new_cols, int new_rows) {
    int pad_width_row = (new_rows - rows) / 2;
    int pad_width_col = (new_cols - cols) / 2;

    // Copy elements from padded back to output within the original bounds of rows and cols
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            output[INDN(i, j, rows)] = padded[INDN(i + pad_width_row, j + pad_width_col, new_rows)];
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

void printmat(double* mat, int size, int stride){
for (int i = 0; i<size; i++){
    for (int j = 0; j<size; j++){
      printf("%lf ", mat[i * (size+stride)+j]);
    }
    printf("\n");
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


int main() {
  int size = next_power_of_2(49);
  
  double inp[25] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};
  double out[size];

  double kernel[9] = {0,0,0,0,1,0,0,0,0};
  double kernelout[size];

  zero_pad(inp, out, 5, 5, 8, 8);
  zero_pad(kernel, kernelout, 3, 3, 8, 8);

  complex cout[size];
  complex ckernelout[size];

  to_complex(out, cout, 8);
  to_complex(kernelout, ckernelout, 8);

  fft2d(cout, -1, 8);
  fft2d(ckernelout, -1, 8);

  complex final[size];

  mmpoint(cout, ckernelout, final, 8);

  fft2d(final, 1, 8);
  fft_shift(final, 8);
  // printcomplexmat(final, 8);

  complex_to_double(final, inp, 8, 5);

  printmat(inp, 5, 0);

  return 0;
}
