// run with ./generate_comp_horner.sh in this directory
// that will compile and run this sample and plot all the generated data
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include "../src/aerobus.h"

using namespace aerobus;  // NOLINT

constexpr size_t NB_POINTS = 400;

template<typename P, typename T, bool compensated>
AEROBUS_DEVICE AEROBUS_INLINE T eval(const T& x) {
    if constexpr (compensated) {
        return P::template compensated_eval<T>(x);
    } else {
        return P::template eval<T>(x);
    }
}

template<typename T>
AEROBUS_DEVICE T exact_large(const T& x) {
    return pow_scalar<T, 5>(0.75 - x) * pow_scalar<T, 11>(1 - x);
}

template<typename T>
AEROBUS_DEVICE T exact_small(const T& x) {
    return pow_scalar<T, 3>(x - 1);
}

// writes x, P(x), exact(x) to plots/file_name and returns the sum of the absolute errors
template<typename P, typename T, bool compensated>
double run(T left, T right, const char *file_name, T (*exact)(const T&)) {
    std::filesystem::create_directories("plots");
    FILE *f = ::fopen((std::string("") + file_name).c_str(), "w+");
    if (f == nullptr) {
        std::perror(file_name);
        return INFINITY;
    }
    double error = 0.0;
    T step = (right - left) / NB_POINTS;
    T x = left;
    for (size_t i = 0; i <= NB_POINTS; ++i) {
        const T value = eval<P, T, compensated>(x);
        ::fprintf(f, "%e %e %e\n", static_cast<double>(x), static_cast<double>(value), static_cast<double>(exact(x)));
        error += std::abs(static_cast<double>(value) - static_cast<double>(exact(x)));
        x += step;
    }
    ::fclose(f);
    return error;
}

int main() {
    bool ok = true;
    {
        // (0.75 - x)^5 * (1 - x)^11
        using P = mul_t<
            pow_t<pq64::val<
                typename q64::template inject_constant_t<-1>,
                q64::val<i64::val<3>, i64::val<4>>>, 5>,
            pow_t<pq64::val<typename q64::template inject_constant_t<-1>, typename q64::one>, 11>
        >;
        using FLOAT = double;
        ok &= run<P, FLOAT, true>(0.68, 1.15, "large_sample_comp_horner.dat", &exact_large)
           <= run<P, FLOAT, false>(0.68, 1.15, "large_sample_horner.dat", &exact_large);
        ok &= run<P, FLOAT, true>(0.74995, 0.75005, "first_root_comp_horner.dat", &exact_large)
           <= run<P, FLOAT, false>(0.74995, 0.75005, "first_root_horner.dat", &exact_large);
        ok &= run<P, FLOAT, true>(0.9935, 1.0065, "second_root_comp_horner.dat", &exact_large)
           <= run<P, FLOAT, false>(0.9935, 1.0065, "second_root_horner.dat", &exact_large);
    }
    {
        // (x - 1) ^ 3: compensated Horner in float against Horner in double
        using P = make_int_polynomial_t<i32, 1, -3, 3, -1>;

        run<P, double, false>(1-0.00005, 1+0.00005, "double.dat", &exact_small);
        run<P, float, true>(1-0.00005, 1+0.00005, "float_comp.dat", &exact_small);
    }
    std::printf(ok ? "compensated Horner is at least as accurate as Horner\n" : "compensated Horner lost accuracy\n");
    return ok ? 0 : 1;
}
