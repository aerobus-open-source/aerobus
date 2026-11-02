// Rational fractions with coefficients in Z/2Z, evaluated in Z/2Z
#include <iostream>
#include "../src/aerobus.h"

using FIELD = aerobus::zpz<2>;
using POLYNOMIALS = aerobus::polynomial<FIELD>;

// x^3 + 2x^2 + 1, with coefficients in Z/2Z: actually x^3 + 1
using P = aerobus::make_int_polynomial_t<FIELD, 1, 2, 0, 1>;
// x^3 + 5x^2 + 8x + 1, with coefficients in Z/2Z: actually x^3 + x^2 + 1
using Q = aerobus::make_int_polynomial_t<FIELD, 1, 5, 8, 1>;

int main() {
    // evaluate in the field itself with value_at_t: P(1) = 1 + 1 = 0 and Q(1) = 1 + 1 + 1 = 1 in Z/2Z
    using one = FIELD::one;
    using p1 = P::value_at_t<one>;
    using q1 = Q::value_at_t<one>;
    // P/Q at 1, computed in Z/2Z
    using f1 = FIELD::div_t<p1, q1>;
    std::cout << "P = " << P::to_string() << ", Q = " << Q::to_string() << std::endl;
    std::cout << "expected (P/Q)(1) = 0" << std::endl;
    std::cout << "computed (P/Q)(1) = " << f1::to_string() << std::endl;
    return f1::v == 0 && q1::v == 1 ? 0 : 1;
}
