#include <iostream>
#include "../src/aerobus.h"

// 32 bits integers as Ring
using RING = aerobus::i32;

// value 2 in this Ring
using X = aerobus::inject_constant_t<RING, 2>;

// Z/2Z as quotient ring
using FIELD = aerobus::Quotient<RING, X>;

int main() {
    // inject 2 in FIELD (Z/2Z as algebraic structure)
    using Y = aerobus::inject_constant_t<FIELD, 2>;

    // get the corresponding numerical value
    constexpr int y = aerobus::get<int, Y>();

    // ensure 2 == 0 in Z/2Z
    std::cout << "expected = 0" << std::endl;
    std::cout << "computed = " << y << std::endl;
    return 0;
}
