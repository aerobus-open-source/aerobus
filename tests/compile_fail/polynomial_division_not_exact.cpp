// must not compile: (x + 1) / (2x) needs 1/2, which is not an integer
#include "aerobus.h"

using P = aerobus::pi64;
using quotient = P::div_t<P::val<aerobus::i64::one, aerobus::i64::one>, P::val<aerobus::i64::val<2>, aerobus::i64::zero>>;

int main() {
    return static_cast<int>(quotient::degree);
}
