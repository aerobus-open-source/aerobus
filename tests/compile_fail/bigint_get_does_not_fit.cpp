// must not compile: 10^19 does not fit in int64_t, get must not return -8446744073709551616 (#51)
#include <cstdint>
#include "aerobus.h"

using A = aerobus::bigint::from_string_t<"10000000000000000000">;

int main() {
    return static_cast<int>(A::get<int64_t>() & 1);
}
