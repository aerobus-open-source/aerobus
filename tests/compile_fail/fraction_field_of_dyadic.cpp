// must not compile: 1/2 would have two forms, the fraction and the value 2^-1
#include "aerobus.h"

using Q = aerobus::FractionField<aerobus::dyadic>;

int main() {
    return static_cast<int>(Q::one::is_zero_v);
}
