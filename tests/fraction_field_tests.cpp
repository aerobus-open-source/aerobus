#include <gtest/gtest.h>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(fraction_field, get) {
    using half = q32::val<i32::one, i32::val<2>>;
    constexpr float x = half::template get<float>();
    EXPECT_EQ(x, 0.5f);

    using thirdhalf = q32::val<i32::val<3>, i32::val<2>>;
    constexpr float y = thirdhalf::template get<float>();
    EXPECT_EQ(y, 1.5f);

    // 3/2 + x / 2
    using polyA = polynomial<q32>::val<half, thirdhalf>;
    constexpr float a = polyA::eval(2.0f);
    EXPECT_EQ(a, 2.5f);

    // 1/2 + x
    using polyB = polynomial<q32>::val<q32::one, half>;
    using F = fpq32::val<polyA, polyB>;
    constexpr float z = F::eval(2.0f);
    EXPECT_EQ(z, 1.0f);
    constexpr float zz = F::eval(-1.0f);
    EXPECT_EQ(zz, -2.0f);
}

TEST(q32, add) {
    {
        // 1/2 + 1/4 = 3/4
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<4>>;

        using c = q32::add_t<a, b>;

        EXPECT_EQ(c::x::v, 3);
        EXPECT_EQ(c::y::v, 4);
    }
    {
        // 1/2 + 1/3 = 5/6
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<3>>;

        using c = q32::add_t<a, b>;

        EXPECT_EQ(c::x::v, 5);
        EXPECT_EQ(c::y::v, 6);
    }
    {
        // -1/2 + 1/2 = 0
        using a = q32::val<i32::val<-1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<2>>;

        using c = q32::add_t<a, b>;
        EXPECT_TRUE(c::is_zero_v);
    }
}

TEST(q32, is_zero) {
    using a = q32::zero;
    EXPECT_TRUE(a::is_zero_v);
    using b = q32::one;
    EXPECT_FALSE(b::is_zero_v);
}

TEST(q32, gt) {
    {
        using a = q32::inject_constant_t<1>;
        using b = q32::zero;
        EXPECT_TRUE((q32::gt_v<a, b>));
    }
    {
        using a = q32::zero;
        using b = q32::inject_constant_t<2>;
        EXPECT_FALSE((q32::gt_v<a, b>));
    }
    {
        using a = q32::zero;
        using b = q32::zero;
        EXPECT_FALSE((q32::gt_v<a, b>));
    }
}

TEST(q32, sub) {
    {
        // 1/2 - 1/4 = 1/4
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<4>>;

        using c = q32::sub_t<a, b>;
        EXPECT_EQ(c::x::v, 1);
        EXPECT_EQ(c::y::v, 4);
    }
    {
        // 1/2 - 1/3 = 1/6
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<3>>;

        using c = q32::sub_t<a, b>;

        EXPECT_EQ(c::x::v, 1);
        EXPECT_EQ(c::y::v, 6);
    }
    {
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<2>>;

        using c = q32::sub_t<a, b>;
        EXPECT_TRUE(c::is_zero_v);
    }
}

TEST(q32, mul) {
    {
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<4>>;

        using c = q32::mul_t<a, b>;

        EXPECT_EQ(c::x::v, 1);
        EXPECT_EQ(c::y::v, 8);
    }
    {
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<3>>;

        using c = q32::mul_t<a, b>;

        EXPECT_EQ(c::x::v, 1);
        EXPECT_EQ(c::y::v, 6);
    }
    {
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<0>, i32::val<2>>;

        using c = q32::mul_t<a, b>;

        EXPECT_TRUE(c::is_zero_v);
    }
}

TEST(q32, div) {
    {
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<4>>;

        using c = q32::div_t<a, b>;

        EXPECT_EQ(c::x::v, 2);
        EXPECT_EQ(c::y::v, 1);
    }
    {
        using a = q32::val<i32::val<1>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<3>>;

        using c = q32::div_t<a, b>;

        EXPECT_EQ(c::x::v, 3);
        EXPECT_EQ(c::y::v, 2);
    }
    {
        using a = q32::val<i32::val<0>, i32::val<2>>;
        using b = q32::val<i32::val<1>, i32::val<2>>;

        using c = q32::div_t<a, b>;

        EXPECT_TRUE(c::is_zero_v);
    }
    {
        using a = q32::val<i32::val<1>, i32::val<1>>;
        using b = q32::val<i32::val<2>, i32::val<1>>;

        using c = q32::div_t<a, b>;


        EXPECT_EQ(c::x::v, 1);
        EXPECT_EQ(c::y::v, 2);
    }
}

TEST(q32, simplify) {
    using A = q32::val<i32::val<2>, i32::val<2>>;
    using B = q32::val<i32::val<1>, i32::val<1>>;
    using C = q32::val<i32::val<-1>, i32::val<-1>>;
    using D = q32::val<i32::val<1>, i32::val<-2>>;
    EXPECT_TRUE((q32::eq_v<q32::simplify_t<A>, q32::one>));
    EXPECT_TRUE((q32::eq_v<q32::simplify_t<B>, q32::one>));
    EXPECT_TRUE((q32::eq_v<q32::simplify_t<C>, q32::one>));
    EXPECT_TRUE((q32::eq_v<q32::simplify_t<D>, q32::val<i32::val<-1>, i32::val<2>>>));
    EXPECT_TRUE((q32::eq_v<q32::simplify_t<q32::sub_t<q32::zero, q32::zero>>, q32::zero>));
}

TEST(q32, eq) {
    using A = q32::val<i32::val<2>, i32::val<2>>;
    using B = q32::val<i32::val<1>, i32::val<1>>;
    using C = q32::val<i32::val<-1>, i32::val<-1>>;
    EXPECT_TRUE((q32::eq_v<A, B>));
    EXPECT_TRUE((q32::eq_v<A, C>));
}

TEST(fraction_field, fraction_field_of_fraction_field_is_same) {
    using qq32 = FractionField<q32>;
    EXPECT_TRUE((std::is_same_v<q32, qq32>));
}

TEST(q32, order_between_positive_fractions) {
    using two_fifths = make_q32_t<2, 5>;
    using half = make_q32_t<1, 2>;
    EXPECT_FALSE((q32::gt_v<two_fifths, half>));
    EXPECT_TRUE((q32::gt_v<half, two_fifths>));
    EXPECT_TRUE((q32::lt_v<two_fifths, half>));
    EXPECT_FALSE((q32::lt_v<half, two_fifths>));
}

TEST(q32, order_with_negative_fractions) {
    using minus_half = make_q32_t<-1, 2>;
    using minus_third = make_q32_t<-1, 3>;
    EXPECT_TRUE((q32::gt_v<minus_third, minus_half>));
    EXPECT_FALSE((q32::gt_v<minus_half, minus_third>));
    EXPECT_TRUE((q32::gt_v<q32::zero, minus_half>));
    EXPECT_FALSE((q32::gt_v<minus_half, q32::zero>));
    EXPECT_TRUE((q32::gt_v<make_q32_t<1, 3>, minus_half>));
    // unreduced input with negative denominator: -2 / -5 = 2/5
    using two_fifths = q32::val<i32::val<-2>, i32::val<-5>>;
    EXPECT_TRUE((q32::gt_v<make_q32_t<1, 2>, two_fifths>));
    EXPECT_TRUE((q32::lt_v<two_fifths, make_q32_t<1, 2>>));
}

TEST(q32, pos) {
    EXPECT_FALSE((q32::pos_v<q32::zero>));
    EXPECT_FALSE((q32::pos_v<q32::val<i32::zero, i32::val<-1>>>));
    EXPECT_TRUE((q32::pos_v<q32::val<i32::val<-1>, i32::val<-2>>>));
    EXPECT_FALSE((q32::pos_v<q32::val<i32::one, i32::val<-2>>>));
}

TEST(fraction_field, rational_functions_equal_up_to_a_constant_factor) {
    using P = polynomial<q32>;
    using two = make_q32_t<2, 1>;
    // (2x + 2) / (4x + 2) and (x + 1) / (2x + 1)
    using a = fpq32::val<P::val<two, two>, P::val<make_q32_t<4, 1>, two>>;
    using b = fpq32::val<P::val<q32::one, q32::one>, P::val<two, q32::one>>;
    EXPECT_TRUE((fpq32::eq_v<a, b>));
    EXPECT_TRUE((std::is_same_v<fpq32::simplify_t<a>, fpq32::simplify_t<b>>));
    EXPECT_TRUE((std::is_same_v<fpq32::add_t<a, fpq32::zero>, fpq32::add_t<b, fpq32::zero>>));
    // simplified denominators are monic
    EXPECT_TRUE((std::is_same_v<typename fpq32::simplify_t<b>::y, P::val<q32::one, make_q32_t<1, 2>>>));
}

TEST(q32, gcd) {
    // in a field every nonzero value divides every value
    using half = make_q32_t<1, 2>;
    EXPECT_TRUE((std::is_same_v<q32::gcd_t<q32::zero, half>, q32::one>));
    EXPECT_TRUE((std::is_same_v<q32::gcd_t<half, make_q32_t<1, 3>>, q32::one>));
    EXPECT_TRUE((std::is_same_v<q32::gcd_t<q32::zero, q32::zero>, q32::zero>));
}

TEST(fraction_field, integral_domains) {
    EXPECT_TRUE((internal::is_integral_domain_v<i32>));
    EXPECT_TRUE((internal::is_integral_domain_v<zpz<7>>));
    EXPECT_FALSE((internal::is_integral_domain_v<zpz<8>>));
    EXPECT_TRUE((internal::is_integral_domain_v<polynomial<zpz<7>>>));
    EXPECT_FALSE((internal::is_integral_domain_v<polynomial<zpz<8>>>));
    EXPECT_TRUE((IsField<FractionField<polynomial<zpz<7>>>>));
}

TEST(fraction_field, every_value_has_one_type) {
    // written in any form, a fraction is the same type; equality is is_same
    using half = q32::val<i32::val<2>, i32::val<4>>;
    static_assert(std::is_same_v<half, q32::val<i32::val<-1>, i32::val<-2>>>);
    static_assert(half::x::v == 1 && half::y::v == 2);
    static_assert(std::is_same_v<q32::val<i32::zero, i32::val<-5>>, q32::zero>);
    static_assert(std::is_same_v<makefraction_t<i32, i32::val<6>, i32::val<3>>, q32::inject_constant_t<2>>);
    // simplify_t is the identity now
    static_assert(std::is_same_v<q32::simplify_t<half>, half>);
    // rational functions: the common factor x + 1 goes away, the denominator is monic
    using PQ = polynomial<q32>;
    using x_plus_1 = PQ::val<q32::one, q32::one>;
    using two_x_plus_2 = PQ::val<q32::inject_constant_t<2>, q32::inject_constant_t<2>>;
    using F = FractionField<PQ>;
    static_assert(std::is_same_v<F::val<x_plus_1, two_x_plus_2>,
                                 F::val<PQ::val<makefraction_t<i32, i32::one, i32::val<2>>>, PQ::one>>);
    static_assert(F::val<x_plus_1, two_x_plus_2>::y::degree == 0);
}
