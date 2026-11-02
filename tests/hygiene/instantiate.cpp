// Compiled with -Wall -Wextra -Wpedantic -Werror: instantiates the main templates of aerobus.h
#include "aerobus.h"

using namespace aerobus;  // NOLINT

int main() {
    using P = make_int_polynomial_t<i64, 1, -3, 2>;
    using Q = make_frac_polynomial_t<i64, 1, 1>;
    using F3 = zpz<3>;
    using R = polynomial<F3>::val<F3::one, F3::val<2>>;
    using D = pi64::derive_t<P>;
    using Z17 = Quotient<i64, i64::val<17>>;
    double x = P::eval(0.5) + P::compensated_eval(0.5) + Q::eval(0.25) + D::eval(1.0)
        + aerobus::sin<i64, 7>::eval(0.1) + aerobus::tan<i64, 5>::compensated_eval(0.1);
    constexpr auto c = P::coeff_at_t<5>::v + R::coeff_at_t<3>::v + pi64::zero::coeff_at_t<2>::v;
    return static_cast<int>(x) + static_cast<int>(c)
        + Z17::div_t<Z17::one, Z17::inject_constant_t<2>>::type::v
        + static_cast<int>(pi64::mul_t<P, P>::degree);
}
