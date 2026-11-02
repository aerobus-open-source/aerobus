#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>

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

TEST(polynomials, eval_in_integer_types) {
    // x^2
    using P = make_int_polynomial_t<i32, 1, 0, 0>;
    EXPECT_EQ((P::eval<int32_t>(1000)), 1000000);
    EXPECT_EQ((P::eval<int64_t>(100000)), 10000000000LL);
    EXPECT_EQ((P::eval<int16_t>(100)), 10000);
    EXPECT_EQ((P::eval<uint64_t>(3)), 9u);
    // H4(10) = 16 * 10^4 - 48 * 10^2 + 12
    EXPECT_EQ((known_polynomials::hermite_phys<4>::eval<int32_t>(10)), 155212);
    EXPECT_TRUE((std::is_same_v<decltype(internal::fma_helper<int32_t>::eval(1, 1, 1)), int32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(internal::fma_helper<int16_t>::eval(1, 1, 1)), int16_t>));
}

TEST(compensated_horner, two_prod_is_error_free_float) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(-1000.0f, 1000.0f);
    int failures = 0;
    for (int i = 0; i < 100000; ++i) {
        const float a = dist(rng);
        const float b = dist(rng);
        float x, y;
        internal::two_prod(a, b, &x, &y);
        // a * b is exact in double, and so is x + y when the transformation is error free
        if (static_cast<double>(x) + static_cast<double>(y) != static_cast<double>(a) * static_cast<double>(b)) {
            ++failures;
        }
    }
    EXPECT_EQ(failures, 0);
}

TEST(compensated_horner, two_prod_is_error_free_double) {
    if constexpr (std::numeric_limits<long double>::digits < 64) {
        GTEST_SKIP() << "needs a long double with at least 64 significant bits";
    } else {
        std::mt19937_64 rng(7);
        std::uniform_real_distribution<double> dist(-1000.0, 1000.0);
        std::uniform_int_distribution<int> small(1, 1023);
        int failures = 0;
        for (int i = 0; i < 100000; ++i) {
            const double a = dist(rng);
            // a has 53 significant bits and b at most 10: their product is exact in long double
            const double b = std::ldexp(static_cast<double>(small(rng)), small(rng) % 16 - 8);
            double x, y;
            internal::two_prod(a, b, &x, &y);
            const long double exact = static_cast<long double>(a) * static_cast<long double>(b);
            if (static_cast<long double>(x) + static_cast<long double>(y) != exact) {
                ++failures;
            }
        }
        EXPECT_EQ(failures, 0);
    }
}

TEST(compensated_horner, more_accurate_than_horner_near_a_multiple_root) {
    // (x - 1)^5, expanded: Horner in float loses every digit near x = 1
    using P = make_int_polynomial_t<i32, 1, -5, 10, -10, 5, -1>;
    std::mt19937 rng(3);
    std::uniform_real_distribution<float> dist(1.05f, 1.1f);
    double worst_compensated = 0.0;
    double worst_horner = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const float x = dist(rng);
        const double exact = std::pow(static_cast<double>(x) - 1.0, 5);
        const double compensated = P::compensated_eval(x);
        const double horner = P::eval(x);
        worst_compensated = std::max(worst_compensated, std::abs(compensated - exact) / std::abs(exact));
        worst_horner = std::max(worst_horner, std::abs(horner - exact) / std::abs(exact));
    }
    EXPECT_LT(worst_compensated, 1e-6) << "Horner relative error: " << worst_horner;
    EXPECT_GT(worst_horner, 100 * worst_compensated);
}

TEST(polynomials, eval_at_infinity) {
    constexpr double inf = std::numeric_limits<double>::infinity();
    // x + 1
    using P = make_int_polynomial_t<i64, 1, 1>;
    EXPECT_EQ(P::eval(inf), inf);
    EXPECT_EQ(P::eval(-inf), -inf);
    EXPECT_EQ(P::eval(std::numeric_limits<float>::infinity()), std::numeric_limits<float>::infinity());
    // -x^2
    using Q = make_int_polynomial_t<i64, -1, 0, 0>;
    EXPECT_EQ(Q::eval(inf), -inf);
    EXPECT_EQ(Q::eval(-inf), -inf);
    // the compensated scheme returns the Horner result when it overflows
    EXPECT_EQ(P::compensated_eval(inf), inf);
    EXPECT_EQ(Q::compensated_eval(-inf), -inf);
}
