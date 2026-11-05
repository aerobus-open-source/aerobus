// must not compile: division by zero
#include "aerobus.h"

using Q = aerobus::dyadic::div_t<aerobus::dyadic::inject_constant_t<7>, aerobus::dyadic::zero>;

int main() {
    return static_cast<int>(Q::exponent);
}
