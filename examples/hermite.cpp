/**
 * Example demonstrating how to use aerobus "pow" function and "hermite" polynomials
 */

#include <cmath>
#include <iostream>
#include "../src/aerobus.h"

namespace standardlib {
    double H3(double x) {
        return 8 * std::pow(x, 3) - 12 * x;
    }

    double H4(double x) {
        return 16 * std::pow(x, 4) - 48 * x * x + 12;
    }
}

namespace aerobuslib {
    double H3(double x) {
        return 8 * aerobus::pow_scalar<double, 3>(x) - 12 * x;
    }

    double H4(double x) {
        return 16 * aerobus::pow_scalar<double, 4>(x) - 48 * x * x + 12;
    }
}

int main() {
    const double h3 = aerobus::known_polynomials::hermite_phys<3>::eval(10.0);
    const double h4 = aerobus::known_polynomials::hermite_phys<4>::eval(10.0);
    const int h4_int = aerobus::known_polynomials::hermite_phys<4>::eval<int>(10);
    std::cout << std::hermite(3, 10) << '=' << standardlib::H3(10) << '=' << h3 << '=' << aerobuslib::H3(10) << '\n'
              << std::hermite(4, 10) << '=' << standardlib::H4(10) << '=' << h4 << '=' << aerobuslib::H4(10)
              << '=' << h4_int << '\n';
    return h3 == std::hermite(3, 10) && h4 == std::hermite(4, 10) && h4_int == 155212 ? 0 : 1;
}
