#include <math.h>
#include <stdio.h>

float quick_dirty_expf(float x);

int main(void) {
    float start = -10.0f, stop=5.0f;
    int steps = 30;
    int i;
    for (i = 0; i < steps; ++i) {
        float x = start + (stop - start) * (float) i / steps;
        float exp_x = quick_dirty_expf(x);
        double golden = exp(x);
        double rel_error = fabs(exp_x - golden) / golden;
        printf("%.5e %.10e %.10e %.4a %.4a\n", x, exp_x, golden, fabs(exp_x - golden), rel_error);
    }

    return 0;
}