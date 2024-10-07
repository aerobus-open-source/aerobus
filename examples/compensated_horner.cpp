// run with ./generate_comp_horner.sh in this directory
// that will compile and run this sample and plot all the generated data
#include "../src/aerobus.h"

using namespace aerobus;  // NOLINT

constexpr size_t NB_POINTS = 400;

template<typename P, typename T, bool compensated>
DEVICE INLINED T eval(const T& x) {
    if constexpr (compensated) {
        return P::template compensated_eval<T>(x);
    } else {
        return P::template eval<T>(x);
    }
}

template<typename T>
DEVICE INLINED T exact(const T& x) {
    return pow_scalar<T, 5>(0.75 - x) * pow_scalar<T, 11>(1 - x);
}

template<typename P, typename T, bool compensated>
void run(T left, T right, const char *file_name) {
    FILE *f = ::fopen(file_name, "w+");
    T step = (right - left) / NB_POINTS;
    T x = left;
    for (size_t i = 0; i <= NB_POINTS; ++i) {
        ::fprintf(f, "%e %e %e\n", x, eval<P, T, compensated>(x), exact(x));
        x += step;
    }
    ::fclose(f);
}

int main() {
    // (0.75 - x)^5 * (1 - x)^11
    using P = mul_t<
        pow_t<pq64, pq64::val<
            typename q64::template inject_constant_t<-1>,
            q64::val<i64::val<3>, i64::val<4>>>, 5>,
        pow_t<pq64, pq64::val<typename q64::template inject_constant_t<-1>, typename q64::one>, 11>
    >;

    ::printf("polynomial = %s\n", P::to_string().c_str());

    using FLOAT = double;
    run<P, FLOAT, false>(0.68, 1.15, "plots/large_sample_horner.dat");
    run<P, FLOAT, true>(0.68, 1.15, "plots/large_sample_comp_horner.dat");

    run<P, FLOAT, false>(0.74995, 0.75005, "plots/first_root_horner.dat");
    run<P, FLOAT, true>(0.74995, 0.75005, "plots/first_root_comp_horner.dat");

    run<P, FLOAT, false>(0.9935, 1.0065, "plots/second_root_horner.dat");
    run<P, FLOAT, true>(0.9935, 1.0065, "plots/second_root_comp_horner.dat");
}
