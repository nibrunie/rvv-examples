#include <math.h>
#include <stddef.h>
#include <assert.h>

#include <bench_poly_mult_utils.h>


/** bit reversal of @p x assuming it is a @p n - bit value */
int bit_reverse(int x, int n) {
    return -1;
}

static int ring_add(ring_t ring, int lhs, int rhs) {
    return (lhs + rhs) % ring.modulo;
}

static int ring_sub(ring_t ring, int lhs, int rhs) {
    return (lhs - rhs) % ring.modulo;
}

static int ring_mul(ring_t ring, int lhs, int rhs) {
    return (lhs * rhs) % ring.modulo;
}

static int ring_square(ring_t ring, int src) {
    return ring_mul(ring, src, src);
}

static int ring_power(ring_t ring, int lhs, int n) {
    assert(n >= 0);
    if (n == 0) return 1;
    else if (n == 1) return lhs;
    else {
        int half_n = n / 2;
        int root = ring_power(ring, lhs, half_n);
        int square = ring_square(ring, root);
        // not constant time (if any doubt remained)
        if (n % 2 == 0) return square;
        else return ring_mul(ring, square, lhs);
    }
}

/** transform @p src into the NTT domain */
void poly_fast_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring, int rootOfUnity, int degree, int start, int stride) {
    if (degree == 0) {
        dst->coeffs[0] = src.coeffs[start];
        dst->modulo = src.modulo;
        dst->degree = degree;
    } else if (degree == 1) {
        int even  = src.coeffs[start];
        int odd   = src.coeffs[start+stride];
        int twiddle = 1;
        twiddle = ring_mul(ring, twiddle, odd);
        dst->coeffs[0] = ring_add(ring, even, twiddle);
        dst->coeffs[1] = ring_sub(ring, even, twiddle);
    } else {
        assert(degree % 2 == 1);
        // re-using destination array as temporary storage for sub NTT
        ntt_t ntt_even = {.coeffs = dst->coeffs, .coeffSize = (dst->coeffSize / 2), .degree = (degree / 2), .modulo = src.modulo};
        ntt_t ntt_odd = {.coeffs = (dst->coeffs + (dst->coeffSize / 2)), .coeffSize = (dst->coeffSize / 2), .degree = (degree / 2), .modulo = src.modulo};
        // recursion
        int root_square = ring_square(ring, rootOfUnity);
        poly_fast_ntt_transform(&ntt_even, src, ring, root_square, degree / 2, start, 2 * stride);
        poly_fast_ntt_transform(&ntt_odd, src, ring, root_square, degree / 2, start + stride, 2 * stride);
        // reconstruction
        int i;
        for (i = 0; i <= degree / 2; ++i) {
            int even  = ntt_even.coeffs[i];
            int odd   = ntt_odd.coeffs[i];
            int twiddle = ring_power(ring, rootOfUnity, i);
            twiddle = ring_mul(ring, twiddle, odd);
            dst->coeffs[i] = ring_add(ring, even, twiddle);
            dst->coeffs[i + (degree + 1) / 2] = ring_sub(ring, even, twiddle);
        }
    }
}

/** transform @p src into the NTT domain */
void poly_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring) {
    assert(dst->degree >= src.degree);

    int d;
    for (d = 0; d <= src.degree; d++) {
        int coeff = 0;
        int k;
        for (k = 0; k <= src.degree; ++k) coeff += ring_mul(ring, src.coeffs[k], ring_power(ring, ring.rootOfUnity, k*d));
        coeff %= ring.modulo;
        dst->coeffs[d] = coeff;
    }
}

void ntt_mul(ntt_t* dst, ntt_t lhs, ntt_t rhs) {
    assert(dst->degree >= lhs.degree && dst->degree >= rhs.degree);

    int d;
    for (d = 0; d <= dst->degree; d++) {
        dst->coeffs[d] = (lhs.coeffs[d] * rhs.coeffs[d]) % dst->modulo;
    }
}

/** transform @p src into the polynomial domain back from NTT domain */
void poly_ntt_inv_transform(polynomial_t* dst, ntt_t src, ring_t ring) {
    assert(dst->degree >= src.degree);

    int d;
    for (d = 0; d <= src.degree; d++) {
        int coeff = 0;
        int k;
        for (k = 0; k <= src.degree; ++k) coeff += ring_mul(ring, src.coeffs[k], ring_power(ring, ring.invRootOfUnity, k*d));
        coeff %= ring.modulo;
        coeff *= ring.invDegree;
        coeff %= ring.modulo;
        dst->coeffs[d] = coeff;
    }
}

/** transform @p src into the polynomial domain back from NTT domain (fast version )*/
void poly_fast_inv_ntt_tranform(polynomial_t* dst, ntt_t src, ring_t ring) {
    assert(dst->degree >= src.degree);

    // inverse NTT transform can be performed by doing NTT with 1/root as initial root of unity
    // and then dividing by the number of coefficients
    poly_fast_ntt_transform(dst, src, ring, ring.invRootOfUnity, src.degree, 0, 1);
    int d;
    for (d = 0; d <= src.degree; d++) {
        dst->coeffs[d] *= ring.invDegree;
        dst->coeffs[d] %= ring.modulo;
    }
}

void poly_mult_ntt(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};
    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    ntt_t ntt_rhs = allocate_poly(rhs.degree, 3329);
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    poly_ntt_transform(&ntt_lhs, lhs, ring);
    poly_ntt_transform(&ntt_rhs, rhs, ring);

    ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_rhs);
    poly_ntt_inv_transform(dst, ntt_lhs_times_rhs, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be free (or statically allocated)
}

void poly_mult_fast_ntt(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};
    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    ntt_t ntt_rhs = allocate_poly(rhs.degree, 3329);
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    poly_fast_ntt_transform(&ntt_lhs, lhs, ring, ring.rootOfUnity, lhs.degree, 0, 1);
    poly_fast_ntt_transform(&ntt_rhs, rhs, ring, ring.rootOfUnity, rhs.degree, 0, 1);

    ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_rhs);

    poly_fast_inv_ntt_tranform(dst, ntt_lhs_times_rhs, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be free (or statically allocated)
}