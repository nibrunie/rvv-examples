#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
  int degree;
  uint32_t modulo;
  uint32_t* coeffs;
} polynomial_t;

static polynomial_t allocate_poly(int degree, int modulo) {
  polynomial_t res;
  res.degree = degree;
  res.modulo = modulo;
  res.coeffs = (uint32_t*) malloc(sizeof(uint32_t) * (degree+1));
  return res;
}

static void randomize_poly(polynomial_t dst) {
  int k;
  for (k = 0; k <= dst.degree; ++k) dst.coeffs[k] = rand() % dst.modulo;
}

static int compare_poly(polynomial_t lhs, polynomial_t rhs) {
  if (lhs.degree != rhs.degree) return -1;
  int k;
  for (k = 0; k < lhs.degree; ++k) {
    if (lhs.coeffs[k] != rhs.coeffs[k]) return -(k+1); // non-zero error code
  };
  return 0; // success
}

/** generic type for a polynomial multiplication implementation */
typedef void(poly_mult_func_t)(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs);

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
    int errors;
} poly_mult_bench_result_t;

static poly_mult_bench_result_t accumulate_bench_result(poly_mult_bench_result_t res, poly_mult_bench_result_t new_result) {
  res.perf_count += new_result.perf_count;

  return res;
}