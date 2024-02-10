#include <math.h>
#include <stdio.h>

float quick_dirty_expf(float x);

int main(void) {
    float start = -10.0f, stop=5.0f;
    int steps = 1000000;
    int i;
    double max_error = 0.;
    double max_error_input = 0.;
    for (i = 0; i < steps; ++i) {
        float x = start + (stop - start) * (float) i / steps;
        float exp_x = quick_dirty_expf(x);
        double golden = exp(x);
        double rel_error = fabs(exp_x - golden) / golden;
        if (rel_error > max_error) {
            max_error = rel_error;
            max_error_input = x;
        }
        // printf("%.5e %.10e %.10e %.4a %.4a\n", x, exp_x, golden, fabs(exp_x - golden), rel_error);
        //printf("%.10e\n", exp_x);
    }
    printf("max relative error is %.5e, it is reached for quick_dirty_expf(%a)=%a vs exp(%a)=%a\n", max_error, max_error_input, quick_dirty_expf(max_error_input), max_error_input, exp(max_error_input));

    return 0;
}