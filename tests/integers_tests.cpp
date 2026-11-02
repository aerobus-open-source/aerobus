#include <gtest/gtest.h>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(i32, eq) {
    using a = i32::inject_constant_t<2>;
    using b = i32::inject_constant_t<3>;
    using c = i32::val<2>;
    EXPECT_TRUE((i32::eq_v<a, c>));
    EXPECT_FALSE((i32::eq_v<a, b>));
}

TEST(integers, gcd) {
    using A = i32::gcd_t<i32::val<12>, i32::val<6>>;
    EXPECT_EQ(A::v, 6);

    using B = i32::gcd_t<i32::val<15>, i32::val<6>>;
    EXPECT_EQ(B::v, 3);

    using C = i32::gcd_t<i32::val<5>, i32::val<3>>;
    EXPECT_EQ(C::v, 1);
}
