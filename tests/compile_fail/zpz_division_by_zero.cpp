// must not compile: division by zero in Z/7Z
#include "aerobus.h"

using Z7 = aerobus::zpz<7>;
using quotient = Z7::div_t<Z7::one, Z7::zero>;

int main() {
    return quotient::v;
}
