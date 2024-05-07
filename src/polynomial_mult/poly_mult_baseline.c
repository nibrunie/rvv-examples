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

/** Basic polynomial addition implementation
 *  dst <- lhs + rhs
*/
void poly_add_baseline(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs)
{
    int d;
    assert(rhs.modulo == lhs.modulo && dst->modulo == rhs.modulo);
    int max_degree = lhs.degree > rhs.degree ? lhs.degree : rhs.degree;
    assert(max_degree <= dst->degree);

    for (d = 0; d <= dst->degree; ++d) {
        int k;
        dst->coeffs[d] = 0;
        if (d <= rhs.degree) {
            dst->coeffs[d] += rhs.coeffs[d];
        } 
        if (d <= lhs.degree) {
            dst->coeffs[d] += lhs.coeffs[d];
        } 
        dst->coeffs[d] %= dst->modulo;
    }
}

/** Basic polynomial in-place mul-sub
 *  dst <- dst - factor . X^power . multiplicand
*/
void poly_in_place_mul_sub(polynomial_t* dst, polynomial_t multiplicand, int factor, int power)
{
    int d;
    assert(dst->modulo == multiplicand.modulo);
    assert(power >= 0);

    for (d = 0; d <= multiplicand.degree; ++d) {
        int dst_degree = d + power;
        if (multiplicand.coeffs[d]) {
            assert(dst_degree <= dst->degree);
            dst->coeffs[dst_degree] -= factor * multiplicand.coeffs[d];
            dst->coeffs[dst_degree] %= dst->modulo;
        }
    }
}

int poly_max_non_zero_degree(polynomial_t poly) {
    int d;
    int max_non_zero_degree = 0;
    for (d = 0; d <= poly.degree; d++) if (poly.coeffs[d] != 0) max_non_zero_degree = d;
    return max_non_zero_degree;
}

/** Polynomial in-place modulo operation
 *  dst <- (dst % modulo)
*/
void poly_in_place_mod(polynomial_t* dst, polynomial_t modulo) {
    int max_non_zero_src_degree = poly_max_non_zero_degree(*dst);
    int max_non_zero_mod_degree = poly_max_non_zero_degree(modulo);
    assert(modulo.coeffs[max_non_zero_mod_degree] == 1);

    while (max_non_zero_src_degree >= max_non_zero_mod_degree) {
        int factor = dst->coeffs[max_non_zero_src_degree];

        if (dst->coeffs[max_non_zero_src_degree]) {
            // dst <- dst - facto * X^(max_non_zero_src_degree - max_non_zero_mode_degree * modulo
            poly_in_place_mul_sub(dst, modulo, factor, max_non_zero_src_degree - max_non_zero_mod_degree);
            assert(dst->coeffs[max_non_zero_src_degree] == 0);
        }
        max_non_zero_src_degree--;
    }
}


/** Polynomial modulo operation
 *  dst <- (src % modulo)
 */
void poly_mod(polynomial_t* dst, polynomial_t src, polynomial_t modulo)
{
    int d;
    assert(src.modulo == modulo.modulo && dst->modulo == src.modulo);
    int max_non_zero_mod_degree = poly_max_non_zero_degree(modulo);

    assert(dst->degree >= src.degree); // extra space required for intermediary computation
    assert(modulo.coeffs[max_non_zero_mod_degree] == 1); // only highest monomial with value one are supported

    // dst <- src (initial copy)
    for (d = 0; d <= dst->degree; d++) {
        if (d <= src.degree) dst->coeffs[d] = src.coeffs[d];
        else dst->coeffs[d] = 0;
    }

    poly_in_place_mod(dst, modulo);
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
        for (k = 0; k <= (d <= lhs.degree ? d : lhs.degree); ++k) dst->coeffs[d] += (lhs.coeffs[k] * rhs.coeffs[d - k]);
        // factorizing modulo accross the full accumulation (assume accumulator is wide enough to accomodate the un-moduloed
        // sum without wrapping around)
        dst->coeffs[d] %= lhs.modulo;
    }
    for (d = rhs.degree + 1; d <= lhs.degree + rhs.degree; ++d) {
        int k;
        dst->coeffs[d] = 0;
        for (k = d - rhs.degree; k <= (d <= lhs.degree ? d : lhs.degree); ++k) dst->coeffs[d] += (lhs.coeffs[k] * rhs.coeffs[d - k]);
        // factorizing modulo accross the full accumulation (assume accumulator is wide enough to accomodate the un-moduloed
        // sum without wrapping around)
        dst->coeffs[d] %= lhs.modulo;
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