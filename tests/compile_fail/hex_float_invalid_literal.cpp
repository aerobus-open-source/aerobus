// must not compile: the literal has no binary exponent
#include "aerobus.h"

using F = aerobus::from_hex_float_t<"0x1.8">;

int main() {
    return static_cast<int>(F::is_zero_v);
}
