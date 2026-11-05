// must not compile: a fraction cannot have a zero denominator
#include "aerobus.h"

using F = aerobus::q32::val<aerobus::i32::one, aerobus::i32::zero>;

int main() {
    return F::x::v;
}
