# Examples of use of RISC-V Vector Reduction operations

https://github.com/nibrunie/rvv-examples/tree/master/src/reduction

This repository contains a few toy examples of vector workloads using reduction operations.

Those examples are described in a blog post on https://fprox.substack.com.

# Using Docker as development environment

## Building the image

The examples in this directory require a RISC-V toolchain.
This toolchain can be build inside a Docker container using one of the docker files
present in this repository: https://github.com/nibrunie/rvv-examples/tree/master/riscv-toolchain-hub.Dockerfile

```
# from the top directory of rvv-examples repository
docker build -t riscv:rvv-examples-reduction -f riscv-toolchain-hub.Dockerfile . 
```


## Running the image in a container

```
# from the top directory of rvv-examples repository
docker run  -ti --mount type=bind,source="$(pwd)"/,target=/home/app riscv:rvv-examples-reduction
```

# Building and executing reduction benchmarks

From within the docker container previously built, you can run the reduction benchmarks easily:

```
cd src/reduction
make clean; make SPIKE=/home/riscv_srcs/spike-vka/build/spike PK=/home/riscv_srcs/riscv-gnu-toolchain/build-pk64/pk sim CFLAGS=" -O3 -DCOUNT_INSTRET -DVERBOSE"
```

# Author(s)

Nicolas Brunie
