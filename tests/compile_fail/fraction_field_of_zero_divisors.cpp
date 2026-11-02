// must not compile: Z/8Z has zero divisors (2 * 4 = 0), it has no field of fractions
#include "aerobus.h"

using F = aerobus::FractionField<aerobus::zpz<8>>;

int main() {
    return F::one::x::v;
}
