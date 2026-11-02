// Host-side checks of the WITH_CUDA_FP16 code paths: built by nvcc, runs without a GPU.
#define WITH_CUDA_FP16
#include <cmath>
#include <cstdio>
#include <random>
#include "aerobus.h"

static int failures = 0;

static void check(bool condition, const char* what) {
    if (!condition) {
        std::printf("FAILED: %s\n", what);
        ++failures;
    }
}

// x is the binary16 value nearest to expected
static bool nearest_half(__half x, double expected) {
    const double v = __half2float(x);
    const double ulp = std::ldexp(1.0, std::ilogb(expected) - 10);
    return std::abs(v - expected) <= ulp / 2;
}

int main() {
    using h = aerobus::arithmetic_helpers<__half>;
    const double pi = 3.14159265358979323846;
    check(nearest_half(h::one(), 1.0), "one");
    check(nearest_half(h::half(), 0.5), "half");
    check(nearest_half(h::pi(), pi), "pi");
    check(nearest_half(h::pi_2(), pi / 2), "pi_2");
    check(nearest_half(h::pi_4(), pi / 4), "pi_4");
    check(nearest_half(h::two_pi(), 2 * pi), "two_pi");
    check(nearest_half(h::inv_two_pi(), 1 / (2 * pi)), "inv_two_pi");

    check(__half2float(aerobus::meta_libm<__half>::floor(__float2half(2.5f))) == 2.0f, "floor on the host");

    // two_prod: x + y == a * b exactly (the product of two halves is exact in float);
    // magnitudes stay away from zero, where the rounding error underflows
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> dist(0.5f, 30.0f);
    std::bernoulli_distribution sign(0.5);
    int inexact = 0;
    for (int i = 0; i < 10000; ++i) {
        const __half a = __float2half(sign(rng) ? dist(rng) : -dist(rng));
        const __half b = __float2half(sign(rng) ? dist(rng) : -dist(rng));
        __half x, y;
        aerobus::internal::two_prod(a, b, &x, &y);
        if (__half2float(x) + __half2float(y) != __half2float(a) * __half2float(b)) {
            ++inexact;
        }
    }
    check(inexact == 0, "two_prod<__half> is error free");

    using P = aerobus::make_int_polynomial_t<aerobus::i16, 1, 1>;
    check(__half2float(P::eval(__float2half(2.0f))) == 3.0f, "eval<__half>");

    std::printf("%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
