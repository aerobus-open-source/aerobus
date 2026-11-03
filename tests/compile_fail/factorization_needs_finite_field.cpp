// must not compile: factorization_t needs coefficients in a finite field, and Z/8Z is not one
#include "aerobus.h"

using P = aerobus::make_int_polynomial_t<aerobus::zpz<8>, 1, 0, 7>;
using F = aerobus::factorization_t<P>;

int main() {
    return static_cast<int>(F::factors::length);
}
