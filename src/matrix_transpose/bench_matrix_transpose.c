// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>


/** Display the content of a n x matrix on stdout */
void matrix_dump(float *mat, unsigned n)
{
    size_t i, j;
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) printf(" %.3f ", mat[i * n + j]);
        printf("\n");

    }

}

/** transpose of a n x n matrix
 *
 * @param dst address of destination matrix
 * @param src address of source matrix
 * @param n matrix dimensions
 */
void matrix_transpose(float *dst,
                      float *src,
                      size_t n);

void matrix_transpose_intrinsics_4x4(float *dst,
                                     float *src);

void matrix_transpose_intrinsics(float *dst,
                                 float *src,
                                 size_t n);

void matrix_transpose_intrinsics_loads(float *dst,
                                       float *src,
                                       size_t n); 

void matrix_4x4_transpose_segmented_load_intrinsics(float* dst, float* src);

void matrix_4x4_transpose_segmented_store_intrinsics(float* dst, float* src);

unsigned long matrix_4x4_transpose_vrgather_bench (float* dst, float* src);

unsigned long matrix_4x4_transpose_vslide_bench(float* dst, float* src); 

/** return the value of the instret counter
 *
 *  The instret counter counts the number of retired (executed) instructions.
*/
static unsigned long read_instret(void)
{
  unsigned long instret;
  asm volatile ("rdinstret %0" : "=r" (instret));
  return instret;
}

// Defining a default size fot the inputs and output array
// (can be overloaded during compilation with -DARRAY_SIZE=<value>)
#ifndef MATRIX_SIZE
#define MATRIX_SIZE 4
#endif

float src[MATRIX_SIZE * MATRIX_SIZE];
float dst[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst2[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst3[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst4[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst5[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst6[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst7[MATRIX_SIZE * MATRIX_SIZE] = {0.f};
float dst8[MATRIX_SIZE * MATRIX_SIZE] = {0.f};


int main(void) {
    int i;
    unsigned long start, stop;
    // random initialization of the input arrays
    for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        src[i] = rand() / (float) RAND_MAX;
    }

    printf("source matrix:\n");
    matrix_dump(src, MATRIX_SIZE);

    start = read_instret();
    matrix_transpose(dst, src, MATRIX_SIZE);
    stop = read_instret();

    printf("matrix_transpose result:\n");
    matrix_dump(dst, MATRIX_SIZE);

    printf("matrix_transpose used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = read_instret();
    matrix_transpose_intrinsics_4x4(dst2, src);
    stop = read_instret();

    printf("matrix_transpose_intrinsics_4x4 result:\n");
    matrix_dump(dst2, MATRIX_SIZE);

    printf("matrix_transpose_intrinsics_4x4 used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = read_instret();
    matrix_transpose_intrinsics(dst3, src, MATRIX_SIZE);
    stop = read_instret();

    printf("matrix_transpose_intrinsics result:\n");
    matrix_dump(dst3, MATRIX_SIZE);

    printf("matrix_transpose_intrinsics used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = read_instret();
    matrix_transpose_intrinsics_loads(dst4, src, MATRIX_SIZE);
    stop = read_instret();

    printf("matrix_transpose_intrinsics_loads result:\n");
    matrix_dump(dst4, MATRIX_SIZE);

    printf("matrix_transpose_intrinsics_loads used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = read_instret();
    matrix_4x4_transpose_segmented_load_intrinsics(dst5, src);
    stop = read_instret();

    printf("matrix_4x4_transpose_segmented_load_intrinscs result:\n");
    matrix_dump(dst5, MATRIX_SIZE);

    printf("matrix_4x4_transpose_segmented_load_intrinscs used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = read_instret();
    matrix_4x4_transpose_segmented_store_intrinsics(dst6, src);
    stop = read_instret();

    printf("matrix_4x4_transpose_segmented_store_intrinscs result:\n");
    matrix_dump(dst6, MATRIX_SIZE);

    printf("matrix_4x4_transpose_segmented_store_intrinscs used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = 0;
    stop = matrix_4x4_transpose_vrgather_bench(dst7, src);

    printf("matrix_4x4_transpose_segmented_store_intrinscs result:\n");
    matrix_dump(dst7, MATRIX_SIZE);

    printf("matrix_4x4_transpose_vrgather used %d instruction(s) to tranpose %dx%d=%d element(s) (in registers).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    start = 0;
    stop = matrix_4x4_transpose_vslide_bench(dst8, src);

    printf("matrix_4x4_transpose_segmented_store_intrinscs result:\n");
    matrix_dump(dst8, MATRIX_SIZE);

    printf("matrix_4x4_transpose_vslide used %d instruction(s) to tranpose %dx%d=%d element(s) (in registers).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);
    return 0;
}
