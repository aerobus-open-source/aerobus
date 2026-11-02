// must not compile: 4 is not prime, there is no Conway polynomial CP(4, 2)
#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"

using P = aerobus::ConwayPolynomial<4, 2>::type;

int main() {
    return static_cast<int>(P::degree);
}
