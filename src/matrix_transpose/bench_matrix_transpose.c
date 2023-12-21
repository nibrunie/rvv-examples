// file: bench_matrix_transpose.c
#include <stdio.h>
#include <stdlib.h>

/** transpose of a n x n matrix
 *
 * @param dst address of destination matrix
 * @param src address of source matrix
 * @param n matrix dimensions
 */
void matrix_transpose(float *dst,
                      float *src,
                      size_t n);

/** return the value of the instret counter
 *
 *  The instret counter counts the number of retired (executed) instructions.
*/
unsigned long read_instret(void)
{
  unsigned long instret;
  asm volatile ("rdinstret %0" : "=r" (instret));
  return instret;
}

// Defining a default size fot the inputs and output array
// (can be overloaded during compilation with -DARRAY_SIZE=<value>)
#ifndef MATRIX_SIZE
#define MATRIX_SIZE 16
#endif

float src[MATRIX_SIZE * MATRIX_SIZE];
float dst[MATRIX_SIZE * MATRIX_SIZE] = {0.f};


int main(void) {
    int i;
    // random initialization of the input arrays
    for (i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        src[i] = rand() / (float) RAND_MAX;
    }

    unsigned long start, stop;
    start = read_instret();
    matrix_transpose(dst, src, MATRIX_SIZE);
    stop = read_instret();

    printf("matrix_transpose used %d instruction(s) to tranpose %dx%d=%d element(s).\n",
           stop - start, MATRIX_SIZE, MATRIX_SIZE, MATRIX_SIZE * MATRIX_SIZE);

    return 0;
}
