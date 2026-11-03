// must not compile: an i64 value is not a coefficient of polynomial<q64>
#include "aerobus.h"

using P = aerobus::polynomial<aerobus::q64>::val<aerobus::i64::one, aerobus::q64::one>;

int main() {
    return static_cast<int>(P::degree);
}
