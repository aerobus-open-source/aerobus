#include <gtest/gtest.h>
#include <cstddef>

#include "aerobus.h"

// The degree limits are measured: the highest degree whose computation fits in the integer type.
// Each Taylor series must accept that degree and reject the next one through its requires clause, so that
// the rejection is a constraint failure (testable here) instead of an overflow inside the header.

namespace {
namespace kp = aerobus::known_polynomials;

#define AEROBUS_POLYNOMIAL_PROBE(F) \
    template<size_t d, typename I> constexpr bool has_##F = requires { typename kp::F<d, I>; };
#define AEROBUS_TAYLOR_PROBE(F) \
    template<size_t d, typename I> constexpr bool has_##F = requires { typename aerobus::F<I, d>; };

AEROBUS_POLYNOMIAL_PROBE(legendre)
AEROBUS_TAYLOR_PROBE(exp)
AEROBUS_TAYLOR_PROBE(expm1)
AEROBUS_TAYLOR_PROBE(lnp1)
AEROBUS_TAYLOR_PROBE(atan)
AEROBUS_TAYLOR_PROBE(sin)
AEROBUS_TAYLOR_PROBE(sinh)
AEROBUS_TAYLOR_PROBE(cosh)
AEROBUS_TAYLOR_PROBE(cos)
AEROBUS_TAYLOR_PROBE(geometric_sum)
AEROBUS_TAYLOR_PROBE(asin)
AEROBUS_TAYLOR_PROBE(asinh)
AEROBUS_TAYLOR_PROBE(atanh)
AEROBUS_TAYLOR_PROBE(tan)
AEROBUS_TAYLOR_PROBE(tanh)

#define EXPECT_MAX_DEGREE(F, I, max) \
    static_assert(has_##F<max, aerobus::I>, #F "<" #max ", " #I "> must be accepted"); \
    static_assert(!has_##F<max + 1, aerobus::I>, #F "<" #max " + 1, " #I "> must be rejected")

#define EXPECT_MAX_DEGREES(F, d16, d32, d64) \
    EXPECT_MAX_DEGREE(F, i16, d16); EXPECT_MAX_DEGREE(F, i32, d32); EXPECT_MAX_DEGREE(F, i64, d64)
}  // namespace

TEST(degree_limits, taylor_series) {
    EXPECT_MAX_DEGREES(exp, 7, 12, 20);
    EXPECT_MAX_DEGREES(expm1, 7, 12, 20);
    EXPECT_MAX_DEGREES(sin, 8, 12, 20);
    EXPECT_MAX_DEGREES(sinh, 8, 12, 20);
    EXPECT_MAX_DEGREES(cos, 7, 13, 21);
    EXPECT_MAX_DEGREES(cosh, 7, 13, 21);
    EXPECT_MAX_DEGREES(asin, 8, 12, 20);
    EXPECT_MAX_DEGREES(asinh, 8, 12, 20);
    EXPECT_MAX_DEGREES(tan, 6, 12, 20);
    EXPECT_MAX_DEGREES(tanh, 6, 12, 20);
}

TEST(degree_limits, series_without_overflow_are_not_limited) {
    static_assert(has_lnp1<120, aerobus::i16> && has_atan<120, aerobus::i16>);
    static_assert(has_atanh<120, aerobus::i16> && has_geometric_sum<120, aerobus::i16>);
}

TEST(degree_limits, other_integer_types_are_not_limited) {
    // no measured limit for a user-defined integer ring: the alias is not constrained
    static_assert(has_legendre<3, aerobus::integers<int8_t>>);
}
