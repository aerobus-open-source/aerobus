// must not compile: 2^32 does not fit in int32_t, inject_constant_t must not truncate it to 0 (#50)
#include <cstdint>
#include "aerobus.h"

using A = aerobus::i32::inject_constant_t<int64_t{1} << 32>;

int main() {
    return A::v;
}
