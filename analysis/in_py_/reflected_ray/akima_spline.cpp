#include "/home/murali/Documents/rass/data/cpp_lib/alglib-cpp/src/interpolation.h"
#include <vector>

extern "C" {
    void build_akima_spline(const double* x, const double* y, int n, double* result, int result_size) {
        alglib::real_1d_array ax, ay;
        ax.setcontent(n, x);
        ay.setcontent(n, y);

        alglib::spline1dinterpolant s;
        alglib::spline1dbuildakima(ax, ay, s);

        for (int i = 0; i < result_size; i++) {
            result[i] = alglib::spline1dcalc(s, x[i]);
        }
    }
}

