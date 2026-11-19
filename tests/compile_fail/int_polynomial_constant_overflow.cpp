// must not compile: 5000000000 does not fit in int32_t, it must not become 705032704 (#50)
#include "aerobus.h"

using P = aerobus::make_int_polynomial_t<aerobus::i32, 5000000000>;

int main() {
    return static_cast<int>(P::degree);
}
