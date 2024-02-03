#include <stddef.h>
#include <stdio.h>

/** generic type for a binary32/float softmax implementation */
typedef void(softmax_func_t)(float* dst, float* src, size_t n);

/** return the value of selected perf counter
 * 
 * perf counter is selected through a macro:
 * - defining COUNT_INSTRET selects the instret counter
 *    The instret counter counts the number of retired (executed) instructions.
 * - defining COUNT_CYCLE selects cycle count
*/
static unsigned long read_perf_counter(void)
{
  unsigned long counter_value;
#if defined(COUNT_INSTRET)
#define PERF_METRIC "instruction"
  asm volatile ("rdinstret %0" : "=r" (counter_value));
#elif defined(COUNT_CYCLE)
#define PERF_METRIC "cycle"
  asm volatile ("rdcycle %0" : "=r" (counter_value));
#else
  // instret is also the default
#define PERF_METRIC "instruction"
  asm volatile ("rdinstret %0" : "=r" (counter_value));
#endif
  return counter_value;
}

typedef struct {
    unsigned long perf_count;
    double max_abs_error;
    double max_rel_error;
    double error_norm2;
} softmax_bench_result_t;


/** Display the content of a binary32 n-element array */
static void array_dump_fp32(float *array, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        printf(" %.3f ", array[i]);
    }
    printf("\n");
}

/** Display the content of a binary64 n-element array */
static void array_dump_fp64(double *array, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) {
        printf(" %.3f ", array[i]);
    }
    printf("\n");
}