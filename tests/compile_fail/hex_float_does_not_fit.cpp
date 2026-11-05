// must not compile: 2^-40 does not fit in q32
#include "aerobus.h"

using F = aerobus::from_hex_float_t<"0x1p-40", aerobus::i32>;

int main() {
    return F::x::v;
}
