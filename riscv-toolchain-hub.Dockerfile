#
# RISC-V Development Dockerfile
#
# Inspired by https://github.com/sbates130272/docker-riscv
#
# This Dockerfile creates a container full of lots of useful tools for
# RISC-V development.

# Pull base image 
FROM nibrunie/riscv-toolchain:riscv-toolchain-20231225-squashed

# Set the maintainer
MAINTAINER Nicolas Brunie

# Make a working folder and set the necessary environment variables.
ENV RISCV /opt/riscv
ENV NUMJOBS 1
RUN mkdir -p $RISCV

# Make a directory to download sources and build programs
ENV RISCV_SRCS /home/riscv_srcs/
RUN mkdir -p $RISCV_SRCS

# Add the GNU utils bin folder to the path.
ENV PATH $RISCV/bin:$PATH

# building missing pk
WORKDIR $RISCV_SRCS/riscv-gnu-toolchain
RUN make stamps/build-pk64 -j5
RUN make install

WORKDIR $RISCV_SRCS
RUN git clone https://github.com/nibrunieAtSi5/riscv-isa-sim.git -b vector-crypto-additional spike-vka
WORKDIR $RISCV_SRCS/spike-vka
RUN mkdir -p build
WORKDIR $RISCV_SRCS/spike-vka/build
RUN ../configure && make -j5


# Defaulting to the /home/app directory (where we are going to bind the rvv-examples directory)
WORKDIR /home/app