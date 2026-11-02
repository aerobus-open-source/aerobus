// must not compile: division by zero in the rationals
#include "aerobus.h"

using aerobus::q32;
using quotient = q32::div_t<q32::zero, q32::zero>;

int main() {
    return quotient::x::v;
}
