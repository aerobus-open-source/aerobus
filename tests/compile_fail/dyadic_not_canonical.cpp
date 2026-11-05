// must not compile: the mantissa is even
#include "aerobus.h"

using Q = aerobus::dyadic::val<aerobus::bigint::inject_constant_t<6>, -2>;

int main() {
    return static_cast<int>(Q::exponent);
}
