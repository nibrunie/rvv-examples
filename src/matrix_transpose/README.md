# Matrix Transpose

The examples in this directory implement multiple version of a simple square matrix transpose.
Those examples are designed to illustrate multiple way RISC-V Vector Extension (RVV) can
be used to perform such operation.
Some examples support any square matrix while others only work with a `4 x 4` matrix.

# Running example in a docker container

All the examples in this directory can be built and executed using a docker container built
from https://github.com/nibrunie/rvv-examples/blob/main/riscv-toolchain.Dockerfile.
You can follow the indication of https://github.com/nibrunie/rvv-examples/blob/main/README.md for more 
information.

# How to build and execute

Once inside a proper build/execution environment, you can change the current directory to
be `src/matrix_transpose` and execute a simple `make` command to build and run the examples:

```
make clean
make sim_bench_matrix_transpose
```

It is possible to alterate build configuration with the `EXTRA_CFLAGS` environment variable,
for example to count cycles rather than retired instructions you can run:
```
make clean
make sim_bench_matrix_transpose EXTRA_CFLAGS=-DCOUNT_CYCLE
```