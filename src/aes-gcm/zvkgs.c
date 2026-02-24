// SPDX-FileCopyrightText: Copyright (c) 2022 by Rivos Inc.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Vector Carryless Multiply Accumulate over GHASH Galois-Field routine using
// the proposed Zvkgs instructions (vghshs.vs).
//
// This code was developed to validate the design of the Zvkgs extension,
// understand and demonstrate expected usage patterns.
//
// DISCLAIMER OF WARRANTY:
//  This code is not intended for use in real cryptographic applications,
//  has not been reviewed, even less audited by cryptography or security
//  experts, etc.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS
//  OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
//  IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
//  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
//  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include <inttypes.h>

#include "zvkgs.h"


// zvkgs_vghsh_vs
//
// Performs one step of GHASH function as described in NIST GCM publication.
// It is a thin wrapper around the vghsh.vs instruction from the Zvkgs extension.
//
void zvkgs_vghsh_vs(
    uint128_t Y[1],  // a0
    uint128_t X[1],  // a1
    uint128_t H[1]   // a2
) {
    // 4 * 32b = 128b
    // We use LMUL=4 to enable runs with VLEN=32, as a proof of concept.
    // Once VLEN>=128, we can simply use LMUL=1.
    __asm__ volatile (
        "vsetivli x0, 4, e32, m4, ta, ma\n"

        "vle32.v v0, (a0)\n"
        "vle32.v v4, (a1)\n"
        "vle32.v v8, (a2)\n"

        "vghsh.vs v0, v8, v4\n"
        "vse32.v v0, (a0)\n"
        : // no output
        : [a0] "r" (Y), [a1] "r" (X), [a2] "r" (H) // Input operands: X and H
        : "v0", "v4", "v8", "memory" // Clobbered registers
    );

}


