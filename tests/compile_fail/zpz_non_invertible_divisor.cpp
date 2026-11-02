// must not compile: 2 has no inverse modulo 8
#include "aerobus.h"

using Z8 = aerobus::zpz<8>;
using quotient = Z8::div_t<Z8::one, Z8::val<2>>;

int main() {
    return quotient::v;
}
