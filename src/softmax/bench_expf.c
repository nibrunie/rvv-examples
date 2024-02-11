#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

float quick_dirty_expf(float x);

float quick_dirty_vector_expf(float* dst, float* src, float max_x, size_t n);

float quick_dirty_expf_from_vector(float x) {
    float src[4] = {x, x, x, x};
    float dst[4] = {0};
    quick_dirty_vector_expf(dst, src, 0.f, 4);
    return dst[3];
}

typedef float (exp_implementation_t)(float);

typedef struct {
    exp_implementation_t* func;
    double max_rel_error;
    float max_rel_error_input;
    char label[100];
} bench_exp_t;

int main(void) {
    srand(17);
    float start = -2.0f, stop=2.0f;
    int steps = 128;
    int i;
    double max_error[2] = {0.};
    double max_error_input[2] = {0.};

    bench_exp_t benchmarks[] = {
        (bench_exp_t) {.func = expf,             .label="baseline expf"},
        (bench_exp_t) {.func = quick_dirty_expf, .label="quick_dirty_expf"},
        (bench_exp_t) {.func = quick_dirty_expf_from_vector, .label="quick_dirty_vector_expf"},
    }; 

    int benchId;
    for (benchId = 0; benchId < sizeof(benchmarks) / sizeof(bench_exp_t); benchId ++) {
        benchmarks[benchId].max_rel_error = 0.f;
        benchmarks[benchId].max_rel_error_input = 0.f;
    }


    for (i = 0; i < steps; ++i) {
        float x = start + (stop - start) * rand() / (float) RAND_MAX;
        double golden = exp(x);
    
        for (benchId = 0; benchId < sizeof(benchmarks) / sizeof(bench_exp_t); benchId ++) {
            float result = benchmarks[benchId].func(x);
            double rel_error = (result - golden) / golden;
            if (rel_error > benchmarks[benchId].max_rel_error) {
                benchmarks[benchId].max_rel_error = rel_error;
                benchmarks[benchId].max_rel_error_input = x;
            }
            printf("%.5e ", rel_error);
        }
        printf("\n");
    }
    for (benchId = 0; benchId < sizeof(benchmarks) / sizeof(bench_exp_t); benchId ++) {
        printf("%s:             max relative error is %.5e, it is reached for expf(%a)=%a vs exp(%a)=%a\n",
               benchmarks[benchId].label, benchmarks[benchId].max_rel_error,
               benchmarks[benchId].max_rel_error_input, benchmarks[benchId].func(benchmarks[benchId].max_rel_error_input),
               benchmarks[benchId].max_rel_error_input, exp(benchmarks[benchId].max_rel_error_input));
    }

    return 0;
}