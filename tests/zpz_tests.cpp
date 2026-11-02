#include <gtest/gtest.h>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(zpz, basic_assertions) {
    using Z2Z = zpz<2>;
    EXPECT_EQ((Z2Z::template add_t<typename Z2Z::val<1>, typename Z2Z::val<1>>::v), 0);
    EXPECT_TRUE(Z2Z::is_field);

    using Z4Z = zpz<4>;
    EXPECT_EQ((Z4Z::template add_t<typename Z4Z::val<4>, typename Z4Z::val<12>>::v), 0);
    EXPECT_NE((Z4Z::template add_t<typename Z4Z::val<5>, typename Z4Z::val<13>>::v), 0);
    EXPECT_FALSE(Z4Z::is_field);

    // gcd
    EXPECT_EQ((Z2Z::template gcd_t<typename Z2Z::val<1>, typename Z2Z::val<1>>::v), 1);

    using Z5Z = zpz<5>;
    EXPECT_EQ((Z5Z::template gcd_t<typename Z5Z::val<2>, typename Z5Z::val<3>>::v), 1);
    EXPECT_EQ((Z5Z::template gcd_t<typename Z5Z::val<2>, typename Z5Z::val<4>>::v), 2);
}
