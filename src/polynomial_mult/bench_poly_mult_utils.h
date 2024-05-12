#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
  int n;
  int modulo;
  int rootOfUnity;
  int invRootOfUnity;
  int invDegree;
} ring_t;

/** Structure to held a polynomial
 *  coeffiSize indicates the size of the coefficient array
 *  degree current polynomial degree (must verify (degree + 1) <= coeffSize) 
*/
typedef struct {
  int degree; // degree of the polynomial currently stored
  int modulo;
  int* coeffs;
  size_t coeffSize; // size of coeffs (number of elements)
} polynomial_t;

typedef polynomial_t ntt_t;

static polynomial_t allocate_poly(int degree, int modulo) {
  polynomial_t res;
  res.degree = degree;
  res.coeffSize = degree + 1;
  res.modulo = modulo;
  res.coeffs = (int*) malloc(sizeof(int) * (res.coeffSize));
  return res;
}

static void randomize_poly(polynomial_t dst) {
  int k;
  for (k = 0; k <= dst.degree; ++k) dst.coeffs[k] = rand() % dst.modulo;
}

static int compare_poly(polynomial_t lhs, polynomial_t rhs) {
  if (lhs.degree != rhs.degree) return -1;
  int k;
  for (k = 0; k <= lhs.degree; ++k) {
    if (lhs.coeffs[k] != rhs.coeffs[k]) return -(k+1); // non-zero error code
  };
  return 0; // success
}

static void print_poly(char* label, polynomial_t poly) {
  printf("%s (coeffSize=%zu, degree=%d): ", label, poly.coeffSize, poly.degree);
  if (poly.coeffs[0] != 0) printf("%d", poly.coeffs[0]);
  if (poly.coeffs[1] != 0) printf(" + %d.X", poly.coeffs[1]);
  int i;
  for (i = 2; i <= poly.degree; ++i) {
    if (poly.coeffs[i] != 0) printf(" + %d.X^%d", poly.coeffs[i], i);
  }
  printf("\n");
}

/** generic type for a polynomial multiplication implementation */
typedef void(poly_mult_func_t)(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs);

/** generic type for a polynomial multiplication (with modulo reduction) implementation */
typedef void(poly_mult_mod_func_t)(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo);

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
  res.errors     |= new_result.errors;

  return res;
}


/** generic benchmark wrapper for polynomial multiplication implementation with modulo reduction
 * (lhs * rhs) % modulo
 *
 * @param dst destination polynomial
 * @param lhs left hand side input polynomial
 * @param rhs right hand side input polynomial
 * @param modulo ring modulo polynomial
 * @param func function under test
 * @param golden golden value for destination array (must contain n valid elements)
*/
static poly_mult_bench_result_t poly_mult_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* modulo, poly_mult_mod_func_t func, polynomial_t* golden) {
    unsigned long start, stop;
    poly_mult_bench_result_t bench_result;

    start = read_perf_counter();
    func(dst, *lhs, *rhs, *modulo);
    stop = read_perf_counter();
    bench_result.perf_count = stop - start;

    bench_result.errors = compare_poly(*dst, *golden);

#       ifdef VERY_VERBOSE
        printf("errors: %d\n", bench_result.errors);
#       endif
    return bench_result;
}
