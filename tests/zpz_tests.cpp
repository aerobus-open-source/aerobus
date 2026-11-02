#include <gtest/gtest.h>
#include <type_traits>

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

TEST(zpz, representatives_are_reduced) {
    using Z7 = zpz<7>;
    EXPECT_EQ((Z7::val<-1>::v), 6);
    EXPECT_EQ((Z7::val<-8>::v), 6);
    EXPECT_EQ((Z7::val<13>::v), 6);
    EXPECT_EQ(Z7::val<-1>::to_string(), "6");
    EXPECT_EQ((Z7::val<-1>::get<int>()), 6);
    EXPECT_TRUE((Z7::eq_v<Z7::val<-1>, Z7::val<6>>));
    // one type per residue class
    EXPECT_TRUE((std::is_same_v<Z7::val<-1>, Z7::val<6>>));
    EXPECT_TRUE((std::is_same_v<Z7::val<7>, Z7::zero>));
    EXPECT_TRUE((std::is_same_v<Z7::inject_constant_t<-6>, Z7::one>));
}

TEST(zpz, operations_stay_in_range) {
    using Z7 = zpz<7>;
    // 2 - 5 = -3 = 4
    EXPECT_TRUE((std::is_same_v<Z7::sub_t<Z7::val<2>, Z7::val<5>>, Z7::val<4>>));
    EXPECT_TRUE((std::is_same_v<Z7::add_t<Z7::val<-1>, Z7::one>, Z7::zero>));
    EXPECT_TRUE((std::is_same_v<Z7::mul_t<Z7::val<-1>, Z7::val<-1>>, Z7::one>));
    // comparisons use the representatives in [0, p)
    EXPECT_TRUE((Z7::gt_v<Z7::val<-1>, Z7::val<5>>));
    EXPECT_TRUE((Z7::lt_v<Z7::val<5>, Z7::val<-1>>));
    EXPECT_TRUE((Z7::pos_v<Z7::val<-1>>));
}

TEST(zpz, large_modulus_does_not_overflow) {
    // p = 2^31 - 1: sums and products of representatives overflow int32_t
    using Zp = zpz<2147483647>;
    using minus_one = Zp::val<2147483646>;
    EXPECT_TRUE((std::is_same_v<Zp::add_t<minus_one, minus_one>, Zp::val<2147483645>>));
    EXPECT_TRUE((std::is_same_v<Zp::mul_t<minus_one, minus_one>, Zp::one>));
    EXPECT_TRUE((std::is_same_v<Zp::sub_t<Zp::zero, minus_one>, Zp::one>));
}

TEST(zpz, gcd_is_a_zpz_value) {
    using Z5 = zpz<5>;
    using g = Z5::gcd_t<Z5::val<2>, Z5::val<4>>;
    EXPECT_TRUE((std::is_same_v<typename g::enclosing_type, Z5>));
    EXPECT_TRUE((std::is_same_v<Z5::gcd_t<Z5::zero, Z5::zero>, Z5::zero>));
}

TEST(zpz, polynomials_over_odd_prime_field) {
    // in F3[x], x - 2x = -x = 2x
    using F3 = zpz<3>;
    using PF3 = polynomial<F3>;
    using two_x = PF3::val<F3::val<2>, F3::val<0>>;
    using diff = PF3::sub_t<PF3::X, two_x>;
    EXPECT_TRUE((std::is_same_v<diff, two_x>));
    EXPECT_TRUE((PF3::eq_t<diff, two_x>::value));
}

TEST(zpz, is_field_for_odd_composite_moduli) {
    EXPECT_FALSE(zpz<25>::is_field);
    EXPECT_FALSE(zpz<35>::is_field);
    EXPECT_FALSE(zpz<121>::is_field);
    EXPECT_TRUE(zpz<997>::is_field);
}
