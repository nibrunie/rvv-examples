#include <math.h>
#include <stddef.h>
#include <assert.h>

#include <bench_poly_mult_utils.h>


/** Baseline implementation of polynomial multiplication
 *
 * @param dst destination polynomial
 * @param lhs left hand side input polynomial
 * @param rhs right hand side input polynomial
*/
void poly_mult_baseline(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs)
{
    int d;
    assert(rhs.modulo == lhs.modulo);
    dst->degree = lhs.degree + rhs.degree;
    dst->modulo = lhs.modulo;
    for (d = 0; d <= lhs.degree + rhs.degree; ++d) {
        int k;
        dst->coeffs[d] = 0;
        if (d <= rhs.degree) {
            for (k = 0; k <= (d <= lhs.degree ? d : lhs.degree); ++k) dst->coeffs[d] += (lhs.coeffs[k] * rhs.coeffs[d - k]) % lhs.modulo;
        } else {
            for (k = d - rhs.degree; k <= (d <= lhs.degree ? d : lhs.degree); ++k) dst->coeffs[d] += (lhs.coeffs[k] * rhs.coeffs[d - k]) % lhs.modulo;
        }
    }
}

void poly_mult_scalar_opt(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs)
{
    int d;
    assert(rhs.modulo == lhs.modulo);
    dst->degree = lhs.degree + rhs.degree;
    dst->modulo = lhs.modulo;
    for (d = 0; d <= rhs.degree; ++d) {
        int k;
        dst->coeffs[d] = 0;
        for (k = 0; k <= (d <= lhs.degree ? d : lhs.degree); ++k) dst->coeffs[d] += (lhs.coeffs[k] * rhs.coeffs[d - k]) % lhs.modulo;
    }
    for (d = rhs.degree + 1; d <= lhs.degree + rhs.degree; ++d) {
        int k;
        dst->coeffs[d] = 0;
        for (k = d - rhs.degree; k <= (d <= lhs.degree ? d : lhs.degree); ++k) dst->coeffs[d] += (lhs.coeffs[k] * rhs.coeffs[d - k]) % lhs.modulo;
    }
}

/** generic benchmark wrapper for polynomial multiplication implementation
 *
 * @param dst destination polynomial
 * @param lhs left hand side input polynomial
 * @param rhs right hand side input polynomial
 * @param func function under test
 * @param golden golden value for destination array (must contain n valid elements)
*/
poly_mult_bench_result_t poly_mult_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, poly_mult_func_t func, polynomial_t* golden) {
    unsigned long start, stop;
    poly_mult_bench_result_t bench_result;

    start = read_perf_counter();
    func(dst, *lhs, *rhs);
    stop = read_perf_counter();
    bench_result.perf_count = stop - start;

    bench_result.errors = compare_poly(*dst, *golden);

#       ifdef VERY_VERBOSE
        printf("errors: %d\n", bench_result.errors);
#       endif
    return bench_result;
}


poly_mult_bench_result_t poly_mult_baseline_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, poly_mult_baseline, golden);
}

poly_mult_bench_result_t poly_mult_scalar_opt_bench(polynomial_t* dst, polynomial_t* lhs, polynomial_t* rhs, polynomial_t* golden) {
    return poly_mult_bench(dst, lhs, rhs, poly_mult_scalar_opt, golden);
}