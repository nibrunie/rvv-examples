#include <stdio.h>
#include <stdlib.h>

#include "bench_poly_mult_utils.h"


void poly_mult_baseline(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs);
void poly_in_place_mod(polynomial_t* dst, polynomial_t modulo);

void poly_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring); 
void poly_ntt_inv_transform(polynomial_t* dst, ntt_t src, ring_t ring);

void poly_fast_ntt_transform(ntt_t* dst, polynomial_t src, ring_t ring, int rootOfUnity); 
void poly_fast_inv_ntt_tranform(polynomial_t* dst, ntt_t src, ring_t ring);

void ntt_mul(ntt_t* dst, ntt_t lhs, ntt_t rhs);

int main () {
    polynomial_t pa = allocate_poly(3, 3329);
    polynomial_t pb = allocate_poly(3, 3329);
    polynomial_t pc = allocate_poly(6, 3329);
    polynomial_t mod = allocate_poly(4, 3329);

    // 
    pa.coeffs[0] = pa.coeffs[1] = pa.coeffs[2] = pa.coeffs[3] = 1;
    pb.coeffs[0] = pb.coeffs[3] = 1;
    pb.coeffs[1] = pb.coeffs[2] = 0;
    mod.coeffs[4] = 1;
    mod.coeffs[0] = (3329 - 1);
    mod.coeffs[3] = mod.coeffs[2] = mod.coeffs[1] = 0;
    print_poly("pa", pa);
    print_poly("pb", pb);
    print_poly("mod", mod);

    poly_mult_baseline(&pc, pa, pb);
    print_poly("pc", pc);
    poly_in_place_mod(&pc, mod);
    print_poly("pc %% mod", pc);

    printf("NTT example:\n");
    ring_t ring = {.modulo =3329, .invDegree = 2497, .invRootOfUnity = 1729, .rootOfUnity = 1600};
    ntt_t ntt_pa = allocate_poly(3, 3329);
    polynomial_t pa_inv_ntt = allocate_poly(3, 3329);
    poly_ntt_transform(&ntt_pa, pa, ring);
    print_poly("pa's NTT", ntt_pa);
    poly_ntt_inv_transform(&pa_inv_ntt, ntt_pa, ring);
    print_poly("pa's inv NTT", pa_inv_ntt);

    poly_fast_ntt_transform(&ntt_pa, pa, ring, ring.rootOfUnity);
    print_poly("pa's fast  NTT", ntt_pa);
    poly_ntt_inv_transform(&pa_inv_ntt, ntt_pa, ring);
    print_poly("pa's inv NTT (from fast NTT)", pa_inv_ntt);

    ntt_t ntt_pb = allocate_poly(3, 3329);
    // poly_fast_ntt_transform(&ntt_pb, pb, ring, ring.rootOfUnity);
    poly_ntt_transform(&ntt_pb, pb, ring);
    print_poly("pb's NTT", ntt_pb);
    ntt_t ntt_pa_times_pb = allocate_poly(3, 3329); 
    ntt_mul(&ntt_pa_times_pb, ntt_pa, ntt_pb);
    print_poly("(pa * pb)'s NTT", ntt_pa_times_pb);

    polynomial_t inv_ntt_fast_pa_x_pb = allocate_poly(3, 3329);
    poly_fast_inv_ntt_tranform(&inv_ntt_fast_pa_x_pb, ntt_pa_times_pb, ring);
    print_poly("(pa * pb)'s inv NTT (fast)", inv_ntt_fast_pa_x_pb);

    polynomial_t inv_ntt_pa_x_pb = allocate_poly(3, 3329);
    poly_ntt_inv_transform(&inv_ntt_pa_x_pb, ntt_pa_times_pb, ring);
    print_poly("(pa * pb)'s inv NTT", inv_ntt_pa_x_pb);
    return 0;
}