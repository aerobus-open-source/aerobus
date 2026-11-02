#include <gtest/gtest.h>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(concepts, ring) {
    EXPECT_TRUE((aerobus::IsRing<i32>));
    EXPECT_TRUE((aerobus::IsRing<i64>));
    EXPECT_TRUE((aerobus::IsRing<zpz<3>>));
    EXPECT_TRUE((aerobus::IsRing<aerobus::q32>));
    EXPECT_TRUE((aerobus::IsField<aerobus::q32>));
    EXPECT_TRUE((aerobus::IsRing<aerobus::polynomial<i32>>));
    EXPECT_TRUE((aerobus::IsField<aerobus::FractionField<aerobus::polynomial<i32>>>));
}
