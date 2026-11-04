// must not compile: division by zero
#include "aerobus.h"

using Q = aerobus::bigint::div_t<aerobus::bigint::inject_constant_t<7>, aerobus::bigint::zero>;

int main() {
    return static_cast<int>(Q::size);
}
