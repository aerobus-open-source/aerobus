// must not compile: division by zero in Q[i]
#include "aerobus.h"

using Qi = aerobus::Complex<aerobus::q32>;
using quotient = Qi::div_t<Qi::one, Qi::zero>;

int main() {
    return quotient::is_zero_v ? 0 : 1;
}
