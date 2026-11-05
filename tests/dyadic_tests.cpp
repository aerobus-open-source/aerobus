#include <gtest/gtest.h>
#include <cstdint>
#include <string>
#include <type_traits>

#include "aerobus.h"

// dyadic: exact numbers m 2^e. The arithmetic itself is checked against Python in dyadic_generated_tests.cpp.

namespace {
using aerobus::bigint;
using aerobus::dyadic;
using half = dyadic::from_hex_float_t<"0x1p-1">;
using three_quarters = dyadic::from_hex_float_t<"0x3p-2">;
}  // namespace

TEST(dyadic, is_an_euclidean_domain) {
    static_assert(aerobus::IsRing<dyadic>);
    static_assert(aerobus::IsEuclideanDomain<dyadic>);
    static_assert(!aerobus::IsField<dyadic>);
    static_assert(dyadic::zero::is_zero_v && !dyadic::one::is_zero_v);
    static_assert(dyadic::pos_v<dyadic::one> && !dyadic::pos_v<dyadic::zero>);
}

TEST(dyadic, every_value_has_one_type) {
    // 6 2^-2 = 3 2^-1, written by make_t in any form
    static_assert(std::is_same_v<dyadic::make_t<bigint::inject_constant_t<6>, -2>,
                                 dyadic::make_t<bigint::inject_constant_t<3>, -1>>);
    static_assert(std::is_same_v<dyadic::make_t<bigint::inject_constant_t<6>, -2>::mantissa,
                                 bigint::inject_constant_t<3>>);
    static_assert(dyadic::make_t<bigint::inject_constant_t<6>, -2>::exponent == -1);
    static_assert(std::is_same_v<dyadic::make_t<bigint::zero, 40>, dyadic::zero>);
    static_assert(std::is_same_v<dyadic::inject_constant_t<-8>, dyadic::make_t<bigint::inject_constant_t<-1>, 3>>);
}

TEST(dyadic, embeddings) {
    static_assert(std::is_same_v<aerobus::Embed<aerobus::i64, dyadic>::type<aerobus::i64::val<-12>>,
                                 dyadic::inject_constant_t<-12>>);
    static_assert(std::is_same_v<aerobus::Embed<bigint, dyadic>::type<bigint::inject_constant_t<10>>,
                                 dyadic::make_t<bigint::inject_constant_t<5>, 1>>);
}

TEST(dyadic, exact_where_q64_overflows) {
    // a double coefficient: its denominator 2^92 does not fit in int64_t
    using c = dyadic::from_hex_float_t<"0x1.23456789abcdfp-40">;
    static_assert(c::get<double>() == 0x1.23456789abcdfp-40);
    static_assert(c::exponent == -92);
    // exact product, with its 106 significant bits
    using square = dyadic::mul_t<c, c>;
    static_assert(square::exponent == -184);
    static_assert(square::get<double>() == 0x1.23456789abcdfp-40 * 0x1.23456789abcdfp-40);
}

TEST(dyadic, power_of_two_divisions_are_exact) {
    static_assert(std::is_same_v<dyadic::div_t<dyadic::one, dyadic::inject_constant_t<4>>,
                                 dyadic::from_hex_float_t<"0x1p-2">>);
    static_assert(dyadic::mod_t<dyadic::one, dyadic::inject_constant_t<4>>::is_zero_v);
}

TEST(dyadic, polynomial_coefficients) {
    using PD = aerobus::polynomial<dyadic>;
    using P = PD::val<half, three_quarters>;  // x / 2 + 3 / 4
    // (x / 2 + 3 / 4)^2 = x^2 / 4 + 3 x / 4 + 9 / 16
    using S = PD::val<dyadic::from_hex_float_t<"0x1p-2">, three_quarters, dyadic::from_hex_float_t<"0x9p-4">>;
    static_assert(std::is_same_v<PD::mul_t<P, P>, S>);
    static_assert(P::eval(2.0) == 1.75);
    static_assert(S::eval(0.5F) == 1.0F);
    static_assert(std::is_same_v<P::value_at_t<half>, dyadic::from_hex_float_t<"0x1p+0">>);
    // leading coefficient a power of two: the division is exact
    static_assert(std::is_same_v<PD::div_t<S, P>, P>);
    EXPECT_EQ(P::compensated_eval(2.0), 1.75);
}

TEST(dyadic, to_string) {
    EXPECT_EQ(dyadic::zero::to_string(), "0");
    EXPECT_EQ(three_quarters::to_string(), "3 * 2^-2");
    EXPECT_EQ(dyadic::inject_constant_t<-5>::to_string(), "-5");
}
