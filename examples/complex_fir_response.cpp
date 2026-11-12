/**
 * Frequency response of a complex FIR filter, with taps computed exactly by the compiler.
 *
 * The lowpass filter is the binomial smoother H(z) = ((1 + z^-1) / 2)^6, taps 1, 6, 15, 20, 15, 6, 1 over 64:
 * gain 1 at DC, an exact zero at Nyquist. Multiplying tap k by i^k = exp(i pi k / 2) shifts its response by
 * fs / 4: G(z) = H(-i z) is a complex bandpass filter, which passes +fs/4 with gain 1 and rejects -fs/4 exactly.
 * Mixing by fs / 4 is a classic in signal processing because i^k needs no multiplication; here the taps h_k i^k
 * are exact values of Complex<q32>, and the static_asserts below prove both properties of the filter.
 *
 * The frequency response is G at z = exp(i w), a polynomial in u = z^-1 = exp(-i w) evaluated on an array of
 * std::complex<float>. Build with -O3 -march=x86-64-v4 -fcx-limited-range for the AVX-512 loop: without
 * -fcx-limited-range, each complex product checks for NaN and may call __mulsc3, which blocks vectorization.
 * With clang, add -mprefer-vector-width=512 for zmm registers and -ffp-contract=fast for fused multiply-adds.
 */

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <numbers>
#include <type_traits>
#include <utility>
#include "../src/aerobus.h"

namespace {
using Q = aerobus::q32;
using C = aerobus::Complex<Q>;
using PQ = aerobus::polynomial<Q>;
using PC = aerobus::polynomial<C>;

// H as a polynomial in u = z^-1
using half_sum = PQ::val<aerobus::make_q32_t<1, 2>, aerobus::make_q32_t<1, 2>>;  // (1 + u) / 2
using lowpass = aerobus::pow_t<half_sum, 6>;
static_assert(lowpass::degree == 6);
static_assert(std::is_same_v<lowpass::coeff_at_t<0>, aerobus::make_q32_t<1, 64>>);
static_assert(std::is_same_v<lowpass::coeff_at_t<1>, aerobus::make_q32_t<6, 64>>);
static_assert(std::is_same_v<lowpass::coeff_at_t<3>, aerobus::make_q32_t<20, 64>>);
static_assert(std::is_same_v<lowpass::value_at_t<Q::one>, Q::one>);                     // DC: u = 1
static_assert(std::is_same_v<lowpass::value_at_t<Q::inject_constant_t<-1>>, Q::zero>);  // Nyquist: u = -1

// tap k of G is h_k i^k
template<typename P, size_t... k>
auto shift_by_quarter(std::index_sequence<k...>) -> PC::val<
    C::mul_t<C::inject_ring_t<typename P::template coeff_at_t<sizeof...(k) - 1 - k>>,
             aerobus::pow_t<C::i, sizeof...(k) - 1 - k>>...>;

using bandpass = decltype(shift_by_quarter<lowpass>(std::make_index_sequence<lowpass::degree + 1>{}));
static_assert(std::is_same_v<bandpass::coeff_at_t<1>, C::val<Q::zero, aerobus::make_q32_t<6, 64>>>);
static_assert(std::is_same_v<bandpass::coeff_at_t<2>, C::val<aerobus::make_q32_t<-15, 64>, Q::zero>>);

// w = +pi/2 (+fs/4): u = -i, gain exactly 1; w = -pi/2 (-fs/4): u = i, gain exactly 0
using minus_i = C::val<Q::zero, Q::inject_constant_t<-1>>;
static_assert(std::is_same_v<bandpass::value_at_t<minus_i>, C::one>);
static_assert(std::is_same_v<bandpass::value_at_t<C::i>, C::zero>);
}  // namespace

// the loop to look at in the assembly: h[j] = G(u[j]), with u[j] = exp(-i w_j)
//
// The parts are read and written through float pointers, which std::complex allows (its layout is that of
// float[2]). The plain loop, h[j] = bandpass::eval(u[j]), gives the same results, and GCC vectorizes it, but
// clang 21 does not: it loads and stores each std::complex<float> as one <2 x float> vector, and its loop
// vectorizer gives up on such accesses ("loop not vectorized", without a reason). With float accesses, clang
// vectorizes the loop 16 floats wide with AVX-512.
void frequency_response(const std::complex<float>* u, std::complex<float>* h, size_t n) {
    const float* in = reinterpret_cast<const float*>(u);
    float* out = reinterpret_cast<float*>(h);
    for (size_t j = 0; j < n; ++j) {
        const std::complex<float> r = bandpass::eval(std::complex<float>(in[2 * j], in[2 * j + 1]));
        out[2 * j] = r.real();
        out[2 * j + 1] = r.imag();
    }
}

int main() {
    using cf = std::complex<float>;
    using cd = std::complex<double>;
    constexpr double pi = std::numbers::pi;
    constexpr size_t kPoints = 256;
    cf u[kPoints], h[kPoints];
    double w[kPoints];
    for (size_t j = 0; j < kPoints; ++j) {
        w[j] = -pi + 2 * pi * static_cast<double>(j) / kPoints;
        u[j] = cf(std::polar(1.0, -w[j]));
    }
    frequency_response(u, h, kPoints);

    // G(exp(i w)) = H(exp(i (w - pi/2))) = exp(-3 i (w - pi/2)) cos((w - pi/2) / 2)^6
    double largest_error = 0;
    for (size_t j = 0; j < kPoints; ++j) {
        const double v = w[j] - pi / 2;
        const cd exact = std::polar(std::pow(std::cos(v / 2), 6), -3 * v);
        largest_error = std::fmax(largest_error, std::abs(cd(h[j]) - exact));
    }
    std::printf("   w / pi   |G|\n");
    for (size_t j = 0; j < kPoints; j += kPoints / 8) {
        std::printf("  %+5.2f   %.6f\n", w[j] / pi, std::abs(h[j]));
    }
    std::printf("largest error against the closed form: %.3g\n", largest_error);
    return largest_error <= 1e-6 ? 0 : 1;
}
