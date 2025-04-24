// Copyright 2022 Rivos Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef ZVKGS_H_
#define ZVKGS_H_

#include <stdint.h>

typedef __uint128_t uint128_t;

// Y, X, and H point to 128 bits values, 32b aligned if the processor
// does not support unaligned access.
//
//   Y <- (Y xor X) o H
// where 'o' is the Galois Field Multiplication in GF(2^128).
extern void
zvkgs_vghsh_vs(
    uint128_t Y[1],  // a0
    uint128_t X[1],  // a1
    uint128_t H[1]   // a2
);

#endif  // ZVKGS_H_
