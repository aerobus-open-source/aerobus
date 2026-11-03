// must not compile: the coefficients of the Legendre polynomial of degree 22 overflow int64_t
#include "aerobus.h"

using P = aerobus::known_polynomials::legendre<22>;

int main() {
    return static_cast<int>(P::degree);
}
