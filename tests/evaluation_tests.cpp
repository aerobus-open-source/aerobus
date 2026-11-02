#include <gtest/gtest.h>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(polynomials, eval) {
    // 1 + 2x + 3x^2
    using poly = polynomial<i32>::val<i32::val<3>, i32::val<2>, i32::val<1>>;
    constexpr int v = poly::eval(1);
    EXPECT_EQ(v, 6);
    constexpr float vv = poly::eval(1.0f);
    EXPECT_EQ(vv, 6.0f);

    // 1/2 + 3x/2
    using polyf = polynomial<q32>::val<q32::val<i32::val<3>, i32::val<2>>, q32::val<i32::val<1>, i32::val<2>>>;
    constexpr float vvv = polyf::eval(1.0f);
    EXPECT_EQ(vvv, 2.0f);

    constexpr double vvvv = polyf::eval(-1.0);
    EXPECT_EQ(vvvv, -1.0);
}

TEST(polynomials, compensated_eval) {
    // 1 + 2x + 3x^2
    using poly = make_int_polynomial_t<i32, 3, 2, 1>;
    float vv = poly::compensated_eval(1.0f);
    EXPECT_EQ(vv, 6.0f);

    // 1/2 + 3x/2
    using polyf = polynomial<q32>::val<q32::val<i32::val<3>, i32::val<2>>, q32::val<i32::val<1>, i32::val<2>>>;
    float vvv = polyf::compensated_eval(1.0f);
    EXPECT_EQ(vvv, 2.0f);

    double vvvv = polyf::compensated_eval(-1.0);
    EXPECT_EQ(vvvv, -1.0);
}
