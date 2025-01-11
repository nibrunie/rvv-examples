# Polynomial Multiplication

The examples in this directory implement multiple versions of polynomial multiplication.
Those examples are designed to illustrate multiple ways RISC-V Vector Extension (RVV) can
be used to perform such operation.

# Running example in a docker container

All the examples in this directory can be built and executed using a docker container built
from https://github.com/nibrunie/rvv-examples/blob/main/riscv-toolchain.Dockerfile.
You can follow the indication of https://github.com/nibrunie/rvv-examples/blob/main/README.md for more 
information.

# How to build and execute polynomial multiplication benchmarks

Once inside a proper build/execution environment, you can change the current directory to
be `src/polynomial_mult` and execute a simple `make` command to build and run the examples:

```
make clean
make sim_bench_poly_mult EXTRA_CFLAGS="-DVERBOSE -O3 -DNDEBUG"
```

It is possible to modify the build configuration with the `EXTRA_CFLAGS` environment variable,
for example to count cycles rather than retired instructions you can run:
```
make clean
make sim_bench_poly_mult EXTRA_CFLAGS="-DVERBOSE -O3 -DNDEBUG -DCOUNT_CYCLE"
```

You can also reduce the verbosity level by removing `-DVERBOSE` and
you can change the number of tests executed for each benchmark wit `-DNUM_TESTS=<num-of-tests>`.

## How to run polynomial multiplication benchmark in the conditions of the optimization experiments

To reproduce the test condition (runner is spike, not used for actual latency evaluation) of the experiments
done to optimize the polynomial multiplication in assembly you can run the following commands:

```
make clean
make sim_bench_poly_mult EXTRA_CFLAGS="-O3 -DLMUL=4 -DNLMUL=2 -DWLMUL=8 -DE32_MASK=8 -DCOUNT_CYCLE -DUSE_SLIDE_PAIR_SWAP_NO -DNUM_TESTS=500  -DNDEBUG -DUSE_VREM_MODULO_NO"
```

# How to build and execute basic polynomial test

This directory contains a basic test evaluating a few functions performing polynomial arithmetic.

Once inside a proper build/execution environment, you can change the current directory to
be `src/polynomial_mult` and execute a simple `make` command to build and run the examples:

```
make clean
make sim_basic_poly_test 
```
