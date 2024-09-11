/**
 * Example demonstrating how to implement a custom taylor serie, 
 * here function is f(x) = exp(exp(x) - 1) for which we know the development : 
 * f(x) = sum{B_n / n! x^n} where B_n is the nth Bell number
 */
#include <cmath>
#include <iostream>
#include <iomanip>
#include "../src/aerobus.h"


template<typename T, size_t i>
struct my_coeff {
    using type = aerobus::makefraction_t<T, aerobus::bell_t<T, i>, aerobus::factorial_t<T, i>>;
};

template<size_t deg>
using F = aerobus::taylor<aerobus::i64, my_coeff, deg>;

int main() {
    constexpr double x = F<15>::eval(0.1);
    double xx = std::exp(std::exp(0.1) - 1);
    std::cout << std::setprecision(18) << x << " == " << xx << std::endl;
}
