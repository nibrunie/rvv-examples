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

/** transform @p src into the NTT domain
 *
 *  The actual polynomial processed by this call is src.coeffs[start::stride]
 *
 *  @param dst destination array for NTT coefficients
 *  @param src input polynomial (original one)
 *  @param ring
 *  @param degree current degree of polynomial input
 *  @param start first index in polynomial coefficient array
 *  @param stride index stride in polynomial coefficient array
*/
void poly_fast_ntt_transform_helper(ntt_t* dst, int* coeffs, ring_t ring, int degree, int stride, int level, int rootPowers[8][64]) {
    if (degree == 0) {
        dst->coeffs[0] = coeffs[0];
    } else if (degree == 1) {
        int even  = coeffs[0];
        int odd   = coeffs[stride];
        dst->coeffs[0] = even + odd;
        dst->coeffs[1] = even - odd;
    } else {
        // re-using destination array as temporary storage for sub NTT
        ntt_t ntt_even = {.coeffs = dst->coeffs, .coeffSize = (dst->coeffSize / 2), .degree = (degree / 2), .modulo = dst->modulo};
        ntt_t ntt_odd = {.coeffs = (dst->coeffs + (dst->coeffSize / 2)), .coeffSize = (dst->coeffSize / 2), .degree = (degree / 2), .modulo = dst->modulo};
        // recursion
        poly_fast_ntt_transform_helper(&ntt_even, coeffs, ring, degree / 2, 2 * stride, level+1, rootPowers);
        poly_fast_ntt_transform_helper(&ntt_odd, coeffs + stride, ring, degree / 2, 2 * stride, level+1, rootPowers);
        // reconstruction
        int i;
        for (i = 0; i <= degree / 2; ++i) {
            int even  = ntt_even.coeffs[i];
            int odd   = ntt_odd.coeffs[i];
            int twiddle = rootPowers[level][i];
            twiddle = ring_mul(ring, twiddle, odd);
            dst->coeffs[i] = ring_add(ring, even, twiddle);
            dst->coeffs[i + (degree + 1) / 2] = ring_sub(ring, even, twiddle);
        }
    }
}

/** Flag to indicate that ring structures (included tabulated root powers) have been initialized */
int ringInit = 0;

/** Shared tabulated arrays of root powers */
int ringPowers[1][8][64];
int ringInvPowers[1][8][64];

void initRootPowerTable(ring_t ring, int rootPowers[8][64], int rootOfUnity) {
    // building a share table of required powers of the root of unity
    rootPowers[0][0] = 1;
    rootPowers[0][1] = rootOfUnity;
    for (int p = 2; p < (64 >> 0); ++p) {
        rootPowers[0][p] = ring_mul(ring, rootPowers[0][p-1], rootPowers[0][1]);
    }
    for (int d = 1; d < 7; d++) {
        rootPowers[d][0] = 1;
        rootPowers[d][1] = ring_square(ring, rootPowers[d-1][1]);
        for (int p = 2; p < (64 >> d); ++p) {
            rootPowers[d][p] = ring_mul(ring, rootPowers[d][p-1], rootPowers[d][1]);
        }
    }
}

ring_t getRing(int degree) {
    assert(degree == 127); // only value currently supported
    ring_t ring = {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};
    // first time initialization of shared tables
    if (!ringInit) {
        initRootPowerTable(ring, ringPowers[0], ring.rootOfUnity);
        initRootPowerTable(ring, ringInvPowers[0], ring.invRootOfUnity);

        ringInit = 1;
    }; 

    return ring;
} 


/** transform @p src into the NTT domain
 *
 *  The actual polynomial processed by this call is src.coeffs[start::stride]
 *
 *  @param dst destination array for NTT coefficients
 *  @param src input polynomial (original one)
 *  @param ring
 *  @param rootOfUnity current root of unity used (might be a power of the ring's root of unity)
*/
void poly_fast_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring, int rootOfUnity) {
    int rootPowers[8][64];
    initRootPowerTable(ring, rootPowers, rootOfUnity);
    poly_fast_ntt_transform_helper(dst, src.coeffs, ring, src.degree, 1, 0, rootPowers);
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
    poly_fast_ntt_transform(dst, src, ring, ring.invRootOfUnity);
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
    ring_t ring = getRing(dst->degree); //{.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};
    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    ntt_t ntt_rhs = allocate_poly(rhs.degree, 3329);
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    poly_ntt_transform(&ntt_lhs, lhs, ring);
    poly_ntt_transform(&ntt_rhs, rhs, ring);

    ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_rhs);
    poly_ntt_inv_transform(dst, ntt_lhs_times_rhs, ring);

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_rhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}

void poly_mult_fast_ntt(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs, polynomial_t modulo) {
    // FIXME: ring structure should be a function argument
    // X^4 - 1 Ring: ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    // X^128 - 1 Ring:
    ring_t ring = getRing(dst->degree); // {.modulo =3329, .invDegree = 3303, .invRootOfUnity = 2522, .rootOfUnity = 33};

    ntt_t ntt_lhs = allocate_poly(lhs.degree, 3329);
    // used for both right-hand-side and destination NTT
    ntt_t ntt_lhs_times_rhs = allocate_poly(dst->degree, 3329); 

    poly_fast_ntt_transform_helper(&ntt_lhs, lhs.coeffs, ring, lhs.degree, 1, 0, ringPowers[0]);
    poly_fast_ntt_transform_helper(&ntt_lhs_times_rhs, rhs.coeffs, ring, rhs.degree, 1, 0, ringPowers[0]);

    ntt_mul(&ntt_lhs_times_rhs, ntt_lhs, ntt_lhs_times_rhs);

    poly_fast_ntt_transform_helper(dst, ntt_lhs_times_rhs.coeffs, ring, ntt_lhs_times_rhs.degree, 1, 0, ringInvPowers[0]);
    // division by the degree
    int d;
    for (d = 0; d <= dst->degree; d++) {
        dst->coeffs[d] *= ring.invDegree;
        dst->coeffs[d] %= ring.modulo;
    }

    // FIXME: ntt_rhs and ntt_lhs's coeffs array should be statically allocated
    free(ntt_lhs.coeffs);
    free(ntt_lhs_times_rhs.coeffs);
}