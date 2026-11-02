#include <gtest/gtest.h>
#include <cstdlib>
#include <type_traits>

#include "aerobus.h"

// the examples and the README use the whole namespace, next to C headers declaring ::div_t
using namespace aerobus;  // NOLINT

TEST(api, generic_division_does_not_clash_with_c_div_t) {
    using A = make_int_polynomial_t<i64, 1, 0, -1>;
    using B = make_int_polynomial_t<i64, 1, -1>;
    EXPECT_TRUE((std::is_same_v<divide_t<A, B>, make_int_polynomial_t<i64, 1, 1>>));
    EXPECT_TRUE((std::is_same_v<divide_t<i32::val<7>, i32::val<2>>, i32::val<3>>));
}
