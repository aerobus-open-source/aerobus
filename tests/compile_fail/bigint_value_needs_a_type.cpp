// must not compile: bigint has no native value type, so factorial_v needs the type of the result
#include "aerobus.h"

constexpr auto F = aerobus::factorial_v<aerobus::bigint, 20>;

int main() {
    return static_cast<int>(F);
}
