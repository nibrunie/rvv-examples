#include <stddef.h>

typedef struct {
    unsigned long perf_count;
    int errors;
} ubench_result_t;

static ubench_result_t accumulate_ubench_result(ubench_result_t res, ubench_result_t new_result) {
  res.perf_count += new_result.perf_count;
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
        .errors = 0 \
    }; \
}\

/** Build a latency benchmark for a 2-operand instruction
 *
 * Latency is measure by building a chain of dependent instructions
 */
#define BENCH_LAT_2OP_INSN(op) \
ubench_result_t bench_lat_##op(size_t n) { \
    long start = read_perf_counter(); \
    for (int i = 0; i < n / 16; i++) { \
        asm volatile( \
            "li a0, 3\n" \
            "li a1, 0x1337beef\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
            #op " a1, a1, a0\n" \
        : \
        : \
        : "a0", "a1"\
        ); \
    } \
    long stop = read_perf_counter(); \
    return (ubench_result_t){ \
        .perf_count = (stop - start), \
        .errors = 0 \
    }; \
}\


#define BENCH_THROUGHPUT_2OP_INSN(op) \
ubench_result_t bench_throughput_##op(size_t n) { \
    long start = read_perf_counter(); \
    for (int i = 0; i < n / 13; i++) { \
        asm volatile( \
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
        : \
        : \
        : "a0", "a1", "a2", "a3", "a4", \
          "a5", "a6", "a7", \
          "t0", "t1", "t2", "t3", "t4", \
          "t5", "t6" \
        ); \
    } \
    long stop = read_perf_counter(); \
    return (ubench_result_t){ \
        .perf_count = (stop - start), \
        .errors = 0 \
    }; \
}\

#define BENCH_LAT_INSN_TC(op) (ubench_t){.bench = bench_lat_##op,                .label= "latency benchmark:    " #op}
#define BENCH_THROUGHPUT_INSN_TC(op) (ubench_t){.bench = bench_throughput_##op,  .label= "throughput benchmark: " #op}