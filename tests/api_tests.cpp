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

TEST(api, sixteen_bit_integers_without_cuda) {
    EXPECT_EQ((i16::add_t<i16::val<3>, i16::val<4>>::v), 7);
    EXPECT_EQ((make_q16_t<1, 2>::get<float>()), 0.5f);
    using P = make_int_polynomial_t<i16, 1, 1>;
    EXPECT_EQ(P::eval(2.0f), 3.0f);
    EXPECT_TRUE((std::is_same_v<q16, FractionField<i16>>));
}
