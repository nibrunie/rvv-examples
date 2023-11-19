// file: bench_vector_add.c
#include <stdio.h>
#include <stdlib.h>

/** vector addition
 *
 * @param dst address of destination array
 * @param lhs address of left hand side operand array
 * @param rhs address of right hand side operand array
 * @param avl application vector length (array size)
 */
void vector_add(float *dst,
                float *lhs,
                float *rhs,
                size_t avl);

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
#ifndef ARRAY_SIZE
#define ARRAY_SIZE 1024
#endif

float lhs[ARRAY_SIZE];
float rhs[ARRAY_SIZE];
float dst[ARRAY_SIZE] = {0.f};


int main(void) {
    int i;
    // random initialization of the input arrays
    for (i = 0; i < ARRAY_SIZE; ++i) {
        lhs[i] = rand() / (float) RAND_MAX;
        rhs[i] = rand() / (float) RAND_MAX;
    }

    unsigned long start, stop;
    start = read_instret();
    vector_add(dst, lhs, rhs, ARRAY_SIZE);
    stop = read_instret();

    printf("vector_add_intrinsics used %d instruction(s) to evaluate %d element(s).\n", stop - start, ARRAY_SIZE);

    return 0;
}
