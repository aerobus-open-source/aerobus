#include <gtest/gtest.h>
#include <type_traits>
#include <utility>

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

template<int32_t p, int32_t a, int32_t b>
constexpr bool division_round_trips() {
    using Z = zpz<p>;
    using A = typename Z::template val<a>;
    using B = typename Z::template val<b>;
    return std::is_same_v<typename Z::template mul_t<typename Z::template div_t<A, B>, B>, A>;
}

// (a / b) * b == a for every a and every b != 0
template<int32_t p, size_t... I>
constexpr bool all_divisions_round_trip(std::index_sequence<I...>) {
    return (division_round_trips<p, static_cast<int32_t>(I / (p - 1)), static_cast<int32_t>(1 + I % (p - 1))>()
            && ...);
}

TEST(zpz, division_multiplies_by_the_inverse) {
    using Z7 = zpz<7>;
    // 2 * 4 = 8 = 1
    EXPECT_TRUE((std::is_same_v<Z7::div_t<Z7::one, Z7::val<2>>, Z7::val<4>>));
    EXPECT_TRUE((all_divisions_round_trip<7>(std::make_index_sequence<7 * 6>{})));
    EXPECT_TRUE((all_divisions_round_trip<13>(std::make_index_sequence<13 * 12>{})));
    // division is exact: the remainder is zero
    EXPECT_TRUE((std::is_same_v<Z7::mod_t<Z7::val<5>, Z7::val<3>>, Z7::zero>));
}

TEST(zpz, division_by_units_of_a_composite_modulus) {
    using Z8 = zpz<8>;
    // 3 * 3 = 9 = 1 and 7 * 7 = 49 = 1
    EXPECT_TRUE((std::is_same_v<Z8::div_t<Z8::one, Z8::val<3>>, Z8::val<3>>));
    EXPECT_TRUE((std::is_same_v<Z8::div_t<Z8::val<5>, Z8::val<7>>, Z8::val<3>>));
}

TEST(zpz, polynomial_division_over_prime_field) {
    // in F3[x]: x + 1 = 2 * (2x + 1) + 2
    using F3 = zpz<3>;
    using PF3 = polynomial<F3>;
    using A = PF3::val<F3::one, F3::one>;
    using B = PF3::val<F3::val<2>, F3::one>;
    EXPECT_TRUE((std::is_same_v<PF3::div_t<A, B>, PF3::inject_constant_t<2>>));
    EXPECT_TRUE((std::is_same_v<PF3::mod_t<A, B>, PF3::inject_constant_t<2>>));
}

TEST(zpz, is_field_for_odd_composite_moduli) {
    EXPECT_FALSE(zpz<25>::is_field);
    EXPECT_FALSE(zpz<35>::is_field);
    EXPECT_FALSE(zpz<121>::is_field);
    EXPECT_TRUE(zpz<997>::is_field);
}
