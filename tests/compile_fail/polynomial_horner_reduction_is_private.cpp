// must not compile: horner_reduction_t is an implementation detail of value_at_t
#include "aerobus.h"

using R = aerobus::pi64::horner_reduction_t<aerobus::pi64::X>;

int main() {
    return 0;
}
