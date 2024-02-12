# Softmax

The examples in this directory implement multiple version of a softmax activation function on a 1D array.
Those examples are designed to illustrate multiple way RISC-V Vector Extension (RVV) can
be used to perform such operation.

# Running example in a docker container

All the examples in this directory can be built and executed using a docker container built
from https://github.com/nibrunie/rvv-examples/blob/main/riscv-toolchain.Dockerfile.
You can follow the indication of https://github.com/nibrunie/rvv-examples/blob/main/README.md for more 
information.

# How to build and execute softmax benchmarks

Once inside a proper build/execution environment, you can change the current directory to
be `src/softmax` and execute a simple `make` command to build and run the examples:

```
make clean
make sim_bench_softmax EXTRA_CFLAGS="-DVERBOSE"
```

It is possible to modify the build configuration with the `EXTRA_CFLAGS` environment variable,
for example to count cycles rather than retired instructions you can run:
```
make clean
make sim_bench_softmax EXTRA_CFLAGS="-DVERBOSE -DCOUNT_CYCLE"
```

You can also reduce the verbosity level by removing `-DVERBOSE` and you can change the number of tests executed for each benchmark wit `-DNUM_TESTS=<num-of-tests>`.

# How to build and execute exponential  accuracy benchmarks

This directory contains an accuracy test bench for various scalar implementation of the exponential function.

Once inside a proper build/execution environment, you can change the current directory to
be `src/softmax` and execute a simple `make` command to build and run the examples:

```
make clean
make sim_bench_softmax EXTRA_CFLAGS="-DVERBOSE"
```