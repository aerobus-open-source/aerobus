// Built against an installed aerobus by tests/package/check_package.cmake
#include <aerobus.h>

int main() {
    // (x + 1)^2 evaluated at 2
    using P = aerobus::make_int_polynomial_t<aerobus::i64, 1, 2, 1>;
    return P::eval(2) == 9 ? 0 : 1;
}
