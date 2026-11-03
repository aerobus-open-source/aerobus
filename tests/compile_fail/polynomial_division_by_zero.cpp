// must not compile: division by the zero polynomial
#include "aerobus.h"

using P = aerobus::pi64;
using quotient = P::div_t<P::X, P::zero>;

int main() {
    return static_cast<int>(quotient::degree);
}
