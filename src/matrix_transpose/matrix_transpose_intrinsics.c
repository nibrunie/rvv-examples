#include <riscv_vector.h>
#include <stddef.h>

/** transpose of a n x n matrix
 *
 * @param dst address of destination matrix
 * @param src address of source matrix
 * @param n matrix dimensions
 */
void matrix_transpose(float *dst,
                      float *src,
                      size_t n) 
{
    size_t i, j;
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j) dst[i * n + j] = src[j * n + i];
};