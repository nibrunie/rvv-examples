#include <stddef.h>
#include "riscv_vector.h"

typedef struct {
    unsigned long perf_count;
    int elt_per_op;
    int errors;
} ubench_result_t;

static ubench_result_t accumulate_ubench_result(ubench_result_t res, ubench_result_t new_result) {
  res.perf_count += new_result.perf_count;
  res.elt_per_op  = new_result.elt_per_op;
  res.errors     |= new_result.errors;

  return res;
}

/** function type for micro-benchmarks */
typedef ubench_result_t (bench_func_t)(size_t);

/** Descriptor structure for softmax benchmark */
typedef struct {
    bench_func_t* bench;
    ubench_result_t result;
    char label[256];
} ubench_t;

#define BENCH_0OP_INSN(op) \
ubench_result_t bench_lat_##op(size_t n) { \
    long start = read_perf_counter(); \
    for (int i = 0; i < n / 16; i++) { \
        asm volatile( \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
            #op "\n" \
        : \
        : \
        : \
        ); \
    } \
    long stop = read_perf_counter(); \
    return (ubench_result_t){ \
        .perf_count = (stop - start), \
        .elt_per_op = 1, \
        .errors = 0 \
    }; \
}\

/** Build a latency benchmark for a 2-operand instruction
 *
 * Latency is measure by building a chain of dependent instructions
 */
#define BENCH_LAT_2OP_INSN(op) \
ubench_result_t bench_lat_##op(size_t n) { \
    size_t cnt = n / 16; \
    long start = read_perf_counter(); \
    asm volatile( \
        "li t0, 7\n" \
        "li t2, 0xcafebebe1337beef\n" \
        "mv t1, t2\n" \
    "1:\n" \
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
	    "add t1, t2, t1\n"\
        #op " t1, t1, t0\n" \
        #op " t1, t1, t0\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : "t0", "t1", "t2"\
    ); \
    long stop = read_perf_counter(); \
    long delta = stop - start; \
    cnt = n / 16; \
    start = read_perf_counter(); \
    asm volatile( \
        "li t0, 7\n" \
        "li t1, 0xcafebebe1337beef\n" \
        "li t2, 0xcafebebe1337beef\n" \
    "1:\n" \
	    "add t1, t2, t1\n"\
	    "add t1, t2, t1\n"\
	    "add t1, t2, t1\n"\
	    "add t1, t2, t1\n"\
	    "add t1, t2, t1\n"\
	    "add t1, t2, t1\n"\
	    "add t1, t2, t1\n"\
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : "t0", "t1", "t2"\
    ); \
    stop = read_perf_counter(); \
    long correction = (stop - start); \
    return (ubench_result_t){ \
        .perf_count = (delta - correction), \
        .elt_per_op = 1, \
        .errors = 0 \
    }; \
}\


#define BENCH_THROUGHPUT_2OP_INSN(op) \
ubench_result_t bench_throughput_##op(size_t n) { \
    long start = read_perf_counter(); \
    size_t cnt = n / 13; \
    asm volatile( \
        "li a0, 3\n" \
        "li a1, 0xcafebebe1337beef\n" \
    "1:\n" \
        #op " a2, a1, a0\n" \
        #op " a3, a1, a0\n" \
        #op " a4, a1, a0\n" \
        #op " a5, a1, a0\n" \
        #op " a6, a1, a0\n" \
        #op " a7, a1, a0\n" \
        #op " t0, a1, a0\n" \
        #op " t1, a1, a0\n" \
        #op " t2, a1, a0\n" \
        #op " t3, a1, a0\n" \
        #op " t4, a1, a0\n" \
        #op " t5, a1, a0\n" \
        #op " t6, a1, a0\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : "a0", "a1", "a2", "a3", "a4", \
      "a5", "a6", "a7", \
      "t0", "t1", "t2", "t3", "t4", \
      "t5", "t6" \
    ); \
    long stop = read_perf_counter(); \
    return (ubench_result_t){ \
        .perf_count = (stop - start), \
        .elt_per_op = 1, \
        .errors = 0 \
    }; \
}\

#define BENCH_LAT_INSN_TC(op) (ubench_t){.bench = bench_lat_##op,                .label= "latency " #op }
#define BENCH_THROUGHPUT_INSN_TC(op) (ubench_t){.bench = bench_throughput_##op,  .label= "throughput " #op}


#define BENCH_THROUGHPUT_2OP_FPD_INSN(op) \
ubench_result_t bench_throughput_##op(size_t n) { \
    size_t cnt = n / 13; \
    long start = read_perf_counter(); \
    asm volatile( \
        "li t0, 0x3fcdbeef3fcdbeef\n" \
        "li t1, 0x4abdcafe4abdcafe\n" \
        "fmv.d.x fa0, t0\n" \
        "fmv.d.x fa1, t1\n" \
    "1:\n" \
        #op ".d fa2, fa1, fa0\n" \
        #op ".d fa3, fa1, fa0\n" \
        #op ".d fa4, fa1, fa0\n" \
        #op ".d fa5, fa1, fa0\n" \
        #op ".d fa6, fa1, fa0\n" \
        #op ".d fa7, fa1, fa0\n" \
        #op ".d ft0, fa1, fa0\n" \
        #op ".d ft1, fa1, fa0\n" \
        #op ".d ft2, fa1, fa0\n" \
        #op ".d ft3, fa1, fa0\n" \
        #op ".d ft4, fa1, fa0\n" \
        #op ".d ft5, fa1, fa0\n" \
        #op ".d ft6, fa1, fa0\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : "fa0", "fa1", "fa2", "fa3", "fa4", \
      "fa5", "fa6", "fa7", \
      "ft0", "ft1", "ft2", "ft3", "ft4", \
      "ft5", "ft6" \
    ); \
    long stop = read_perf_counter(); \
    cnt = n / 13; \
    long delta = stop - start; \
    start = read_perf_counter(); \
    asm volatile( \
        "li t0, 0x3fcdbeef3fcdbeef\n" \
        "li t1, 0x4abdcafe4abdcafe\n" \
        "fmv.d.x fa0, t0\n" \
        "fmv.d.x fa1, t1\n" \
    "1:\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : "fa0", "fa1", "fa2", "fa3", "fa4", \
      "fa5", "fa6", "fa7", \
      "ft0", "ft1", "ft2", "ft3", "ft4", \
      "ft5", "ft6" \
    ); \
    stop = read_perf_counter(); \
    long correction = stop - start; \
    return (ubench_result_t){ \
        .perf_count = (delta - correction), \
        .elt_per_op = 1, \
        .errors = 0 \
    }; \
}\

/** Build a latency benchmark for a 2-operand floating-point instruction
 *
 * Latency is measure by building a chain of dependent instructions
 */
#define BENCH_LAT_2OP_FPD_INSN(op) \
ubench_result_t bench_lat_##op(size_t n) { \
    size_t cnt = n / 16; \
    long start = read_perf_counter(); \
    asm volatile( \
        "li t0, 0x3fcdbeef3fcdbeef\n" \
        "li t2, 0x4abdcafe4abdcafe\n" \
        "fmv.d.x ft0, t0\n" \
        "fmv.d.x ft1, t2\n" \
        "fmv.d.x ft2, t2\n" \
    "1:\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        #op ".d ft1, ft1, ft0\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : "ft0", "ft1", "ft2"\
    ); \
    long stop = read_perf_counter(); \
    long delta = stop - start; \
    cnt = n / 16; \
    start = read_perf_counter(); \
    asm volatile( \
        "li t0, 0x3fcdbeef3fcdbeef\n" \
        "li t2, 0x4abdcafe4abdcafe\n" \
        "fmv.d.x ft0, t0\n" \
        "fmv.d.x ft1, t2\n" \
        "fmv.d.x ft2, t2\n" \
    "1:\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : \
    ); \
    stop = read_perf_counter(); \
    long correction = (stop - start); \
    return (ubench_result_t){ \
        .perf_count = (delta - correction), \
        .elt_per_op = 1, \
        .errors = 0 \
    }; \
}


#define BENCH_LAT_VEC_INSN_TC(op, LMUL, elt) (ubench_t){.bench = bench_lat_##op##_m##LMUL##_e##elt, .label= "latency " #op " m" #LMUL " e" #elt }

/** Build a latency benchmark for a 2-operand vector instruction
 *
 * Latency is measure by building a chain of dependent instructions
 */
#define BENCH_LAT_2OP_VEC_INSN(op, LMUL, elt) \
ubench_result_t bench_lat_##op##_m##LMUL##_e##elt(size_t n) { \
    size_t cnt = n / 16; \
    size_t vl = 0; \
    long start = read_perf_counter(); \
    asm volatile( \
        "vsetvli t0, x0, e32, m" #LMUL ", ta, ma\n" \
        "li t0, 0x3c003c00\n" \
        "vmv.v.x v8, t0\n" \
        "vsetvli %[vl], x0, e" #elt ", m" #LMUL ", ta, ma\n" \
        "vid.v v16\n" \
        "vor.vv v16, v16, v8\n" \
        "vid.v v24\n" \
        "vor.vv v24, v24, v8\n" \
    "1:\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        #op ".vv v24, v16, v24\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt), [vl]"+r"(vl) \
    : \
    : "t0" \
    ); \
    long stop = read_perf_counter(); \
    long delta = stop - start; \
    cnt = n / 16; \
    start = read_perf_counter(); \
    asm volatile( \
        "vsetvli t0, x0, e" #elt ", m" #LMUL ", ta, ma\n" \
        "vid.v v16\n" \
        "vid.v v24\n" \
    "1:\n" \
        "addi %[cnt], %[cnt], -1\n" \
        "bnez %[cnt], 1b\n" \
    : [cnt]"+r"(cnt) \
    : \
    : \
    ); \
    stop = read_perf_counter(); \
    long correction = stop - start; \
    return (ubench_result_t){ \
        .perf_count = (delta - correction), \
        .elt_per_op = vl, \
        .errors = 0 \
    }; \
}
