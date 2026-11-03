// must not compile: 'x' is not a decimal digit
#include "aerobus.h"

using N = aerobus::bigint::from_string_t<"12x4">;

int main() {
    return static_cast<int>(N::size);
}
