#include <stdio.h>
#include <stdlib.h>

#include "bench_poly_mult_utils.h"


void poly_mult_baseline(polynomial_t* dst, polynomial_t lhs, polynomial_t rhs);
void poly_in_place_mod(polynomial_t* dst, polynomial_t modulo);

int main () {
    polynomial_t pa = allocate_poly(3, 3329);
    polynomial_t pb = allocate_poly(3, 3329);
    polynomial_t pc = allocate_poly(6, 3329);
    polynomial_t mod = allocate_poly(4, 3329);

    // 
    pa.coeffs[0] = pa.coeffs[1] = pa.coeffs[2] = pa.coeffs[3] = 1;
    pb.coeffs[0] = pb.coeffs[3] = 1;
    pb.coeffs[1] = pb.coeffs[2] = 0;
    mod.coeffs[4] = mod.coeffs[0] = 1;
    mod.coeffs[3] = mod.coeffs[2] = mod.coeffs[1] = 0;
    print_poly("pa", pa);
    print_poly("pb", pb);
    print_poly("mod", mod);

    poly_mult_baseline(&pc, pa, pb);
    print_poly("pc", pc);
    poly_in_place_mod(&pc, mod);
    print_poly("pc %% mod", pc);
    return 0;
}