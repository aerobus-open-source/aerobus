#include <gtest/gtest.h>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(polynomials, simplify) {
    using poly1 = polynomial<i32>::val<i32::val<0>, i32::val<1>, i32::val<2>>;
    using simplified1 = polynomial<i32>::simplify_t<poly1>;
    using expected1 = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
    EXPECT_TRUE((std::is_same_v<expected1, simplified1>));

    using poly2 = polynomial<i32>::val<i32::val<12>>;
    using simplified2 = polynomial<i32>::simplify_t<poly2>;
    EXPECT_TRUE((std::is_same_v<poly2, simplified2>));

    using poly3 = polynomial<i32>::val<i32::val<0>, i32::val<0>, i32::val<1>, i32::val<2>>;
    using simplified3 = polynomial<i32>::simplify_t<poly3>;
    using expected3 = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
    EXPECT_TRUE((std::is_same_v<expected3, simplified3>));
}

TEST(polynomials, construction_strips_leading_zeros) {
    using P = polynomial<i32>::val<i32::zero, i32::zero, i32::one, i32::val<2>>;
    EXPECT_EQ(P::degree, 1);
    EXPECT_TRUE((std::is_same_v<P, polynomial<i32>::val<i32::one, i32::val<2>>>));
    EXPECT_TRUE((std::is_same_v<polynomial<i32>::val<i32::zero, i32::zero>, polynomial<i32>::zero>));
    // a zero leading coefficient that is not the type Ring::zero
    using Q = polynomial<q32>::val<q32::val<i32::zero, i32::val<2>>, q32::one>;
    EXPECT_EQ(Q::degree, 0);
    EXPECT_TRUE((std::is_same_v<Q, polynomial<q32>::one>));
}

TEST(polynomials, eq_compares_values) {
    using PQ = polynomial<q32>;
    // x / 2 + 1, once with an unreduced coefficient
    using A = PQ::val<q32::val<i32::val<2>, i32::val<4>>, q32::one>;
    using B = PQ::val<make_q32_t<1, 2>, q32::one>;
    EXPECT_TRUE((PQ::eq_t<A, B>::value));
    EXPECT_FALSE((PQ::eq_t<A, PQ::X>::value));
}

TEST(polynomials, products_are_simplified_over_rings_with_zero_divisors) {
    // in Z/4Z[x]: (2x + 1)^2 = 4x^2 + 4x + 1 = 1
    using Z4 = zpz<4>;
    using P = polynomial<Z4>::val<Z4::val<2>, Z4::one>;
    EXPECT_TRUE((std::is_same_v<polynomial<Z4>::mul_t<P, P>, polynomial<Z4>::one>));
}

TEST(polynomials, order) {
    using two = pi64::val<i64::val<2>>;
    using one = pi64::one;
    EXPECT_TRUE((pi64::gt_t<two, one>::value));
    EXPECT_FALSE((pi64::gt_t<one, two>::value));
    EXPECT_TRUE((pi64::lt_t<one, two>::value));
    // higher degree is greater
    EXPECT_TRUE((pi64::gt_t<pi64::X, two>::value));
    // same degree and leading coefficient: the next coefficients decide
    using x_plus_1 = make_int_polynomial_t<i64, 1, 1>;
    using x_plus_2 = make_int_polynomial_t<i64, 1, 2>;
    EXPECT_TRUE((pi64::lt_t<x_plus_1, x_plus_2>::value));
    EXPECT_TRUE((pi64::gt_t<x_plus_2, x_plus_1>::value));
    EXPECT_FALSE((pi64::gt_t<x_plus_1, x_plus_1>::value));
    EXPECT_FALSE((pi64::lt_t<x_plus_1, x_plus_1>::value));
}

TEST(polynomials, coeff_at) {
    // 1 + 2x + 3x^2
    using poly = polynomial<i32>::val<i32::val<3>, i32::val<2>, i32::val<1>>;
    using at0 = poly::coeff_at_t<0>;
    using at1 = poly::coeff_at_t<1>;
    using at2 = poly::coeff_at_t<2>;
    using at3 = poly::coeff_at_t<3>;

    EXPECT_TRUE((std::is_same_v<at0, i32::val<1>>));
    EXPECT_TRUE((std::is_same_v<at1, i32::val<2>>));
    EXPECT_TRUE((std::is_same_v<at2, i32::val<3>>));
    EXPECT_TRUE((std::is_same_v<at3, i32::val<0>>));
}

template<typename... coeffs>
using IX = polynomial<i32>::val<coeffs...>;
template<int32_t x>
using Int = i32::val<x>;

template<typename P1, typename P2>
using add_ix = polynomial<i32>::add_t< P1, P2>;
template<typename P1, typename P2>
using sub_ix = polynomial<i32>::sub_t<P1, P2>;

TEST(polynomials, add) {
    {
        // 1 + x
        using P1 = IX<Int<1>, Int<1>>;
        // 1 + x + x^2
        using P2 = IX<Int<1>, Int<1>, Int<1>>;
        // 2 + 2x + x^2
        using A = add_t<P1, P2>;

        EXPECT_EQ((A::coeff_at_t<0>::v), 2);
        EXPECT_EQ((A::coeff_at_t<1>::v), 2);
        EXPECT_EQ((A::coeff_at_t<2>::v), 1);
        EXPECT_EQ(A::degree, 2);
    }

    {
        // 1 + x - x^2
        using P1 = make_int_polynomial_t<i32, -1, 1, 1>;
        // 1 + x + x^2
        using P2 = make_int_polynomial_t<i32, 1, 1, 1>;
        // 2 + 2x
        using A = add_t<P1, P2>;

        EXPECT_EQ((A::coeff_at_t<0>::v), 2);
        EXPECT_EQ((A::coeff_at_t<1>::v), 2);
        EXPECT_EQ((A::coeff_at_t<2>::v), 0);
        EXPECT_EQ(A::degree, 1);
    }
}

TEST(polynomials, derive) {
    {
        // 1 + x
        using P1 = IX<Int<1>, Int<1>>;
        using PP = polynomial<i32>::template derive_t<P1>;
        EXPECT_EQ(PP::degree, 0);
        EXPECT_EQ((PP::coeff_at_t<0>::v), 1);
    }
    {
        // 1 + x + 3x^2
        using P1 = IX<Int<3>, Int<1>, Int<1>>;
        using PP = polynomial<i32>::template derive_t<P1>;
        EXPECT_EQ(PP::degree, 1);
        EXPECT_EQ((PP::coeff_at_t<0>::v), 1);
        EXPECT_EQ((PP::coeff_at_t<1>::v), 6);
    }
    {
        // in Z/2Z
        using Z2Z = zpz<2>;
        // 1 + x + 3x^2 -> 1
        using P1 = polynomial<Z2Z>::template val<Z2Z::template val<3>, Z2Z::template val<1>, Z2Z::template val<1>>;
        using PP = polynomial<Z2Z>::template derive_t<P1>;
        EXPECT_EQ(PP::degree, 0);
        EXPECT_EQ((PP::coeff_at_t<0>::v), 1);
    }
    {
        // in Z/2Z
        using Z2Z = zpz<2>;
        // x^3 + x^2 + x + 1 -> 1 + x^2 (because 3 == 1)
        using P1 = polynomial<Z2Z>::template val<
                            Z2Z::template val<1>,
                            Z2Z::template val<1>,
                            Z2Z::template val<1>,
                            Z2Z::template val<1>>;
        using PP = polynomial<Z2Z>::template derive_t<P1>;
        EXPECT_EQ(PP::degree, 2);
        EXPECT_EQ((PP::coeff_at_t<0>::v), 1);
        EXPECT_EQ((PP::coeff_at_t<1>::v), 0);
        EXPECT_EQ((PP::coeff_at_t<2>::v), 1);
    }
}

TEST(polynomials, sub) {
    {
        // 1 + x
        using P1 = IX<Int<1>, Int<1>>;
        // 1 + x + x^2
        using P2 = IX<Int<1>, Int<1>, Int<1>>;
        // -x2
        using A = sub_ix<P1, P2>;

        EXPECT_EQ((A::coeff_at_t<0>::v), 0);
        EXPECT_EQ((A::coeff_at_t<1>::v), 0);
        EXPECT_EQ((A::coeff_at_t<2>::v), -1);
        EXPECT_EQ(A::degree, 2);
    }

    {
        // 1 + x + x^2
        using P1 = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
        // 1 + x + x^2
        using P2 = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
        // 0
        using A = sub_t<P2, P1>;

        EXPECT_EQ((A::coeff_at_t<0>::v), 0);
        EXPECT_EQ(A::degree, 0);
    }
}

TEST(polynomials, eq) {
    {
        using A = polynomial<i32>::val<i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<1>>;
        EXPECT_TRUE((polynomial<i32>::eq_t<A, B>::value));
    }
    {
        using A = polynomial<i32>::val<i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<2>>;
        EXPECT_FALSE((polynomial<i32>::eq_t<A, B>::value));
    }
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        EXPECT_TRUE((polynomial<i32>::eq_t<A, B>::value));
    }
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<1>, i32::val<2>>;
        EXPECT_FALSE((polynomial<i32>::eq_t<A, B>::value));
    }
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<2>>;
        EXPECT_FALSE((polynomial<i32>::eq_t<A, B>::value));
    }
}

TEST(polynomials, mul) {
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<-1>>;
        using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using mul = mul_t<A, B>;
        using expected = polynomial<i32>::val<i32::val<1>, i32::zero, i32::val<-1>>;
        EXPECT_TRUE((std::is_same_v<expected, mul>));
    }
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using mul = mul_t<A, B>;
        using expected = polynomial<i32>::val<i32::val<1>, i32::val<2>, i32::val<1>>;
        EXPECT_TRUE((std::is_same_v<expected, mul>));
    }
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<2>>;
        using mul = mul_t<A, B>;
        using expected = polynomial<i32>::val<i32::val<2>, i32::val<2>>;
        EXPECT_TRUE((std::is_same_v<expected, mul>));
    }
}

TEST(polynomials, monomial) {
    {
        // 2x^3 + 0 + 0 + 0
        using A = polynomial<i32>::monomial_t<i32::val<2>, 3>;
        EXPECT_EQ((A::coeff_at_t<3>::v), 2);
        EXPECT_EQ((A::coeff_at_t<2>::v), 0);
        EXPECT_EQ((A::coeff_at_t<1>::v), 0);
        EXPECT_EQ((A::coeff_at_t<0>::v), 0);
    }
}

TEST(polynomials, div) {
    // divisibility
    {
        // x^2 - 1
        using A = polynomial<i32>::val<i32::val<1>, i32::val<0>, i32::val<-1>>;
        // x - 1
        using B = polynomial<i32>::val<i32::val<1>, i32::val<-1>>;
        // x + 1
        using Q = aerobus::div_t<A, B>;
        using R = polynomial<i32>::mod_t<A, B>;
        EXPECT_TRUE(R::is_zero_v);
        EXPECT_EQ(Q::degree, 1);
        EXPECT_EQ((Q::coeff_at_t<0>::v), 1);
        EXPECT_EQ((Q::coeff_at_t<1>::v), 1);
    }
    // divide by constant
    {
        using A = polynomial<i32>::val<i32::val<2>, i32::val<2>, i32::val<2>>;
        using B = polynomial<i32>::val<i32::val<2>>;
        using C = polynomial<i32>::div_t<A, B>;
        using R = polynomial<i32>::mod_t<A, B>;
        EXPECT_TRUE(R::is_zero_v);
        EXPECT_EQ(C::degree, 2);
        EXPECT_EQ((C::coeff_at_t<0>::v), 1);
        EXPECT_EQ((C::coeff_at_t<1>::v), 1);
        EXPECT_EQ((C::coeff_at_t<2>::v), 1);
    }
    // no divisibility
    {
        using A = polynomial<i32>::val<i32::val<1>, i32::val<1>, i32::val<1>>;
        using B = polynomial<i32>::val<i32::val<1>, i32::val<1>>;
        using C = polynomial<i32>::div_t<A, B>;
        using R = polynomial<i32>::mod_t<A, B>;
        EXPECT_EQ(C::degree, 1);
        EXPECT_EQ((C::coeff_at_t<0>::v), 0);
        EXPECT_EQ((C::coeff_at_t<1>::v), 1);
        EXPECT_EQ(R::degree, 0);
        EXPECT_EQ(R::aN::v, 1);
    }
    // float divisibility
    {
        // x2 -1
        using A = polynomial<q32>::val<q32::one, q32::zero, q32::val<i32::val<-1>, i32::val<1>>>;
        // 2x + 2
        using B = polynomial<q32>::val<q32::val<i32::val<2>, i32::one>, q32::val<i32::val<2>, i32::one>>;
        using Q = polynomial<q32>::div_t<A, B>;

        EXPECT_EQ(Q::degree, 1);
        EXPECT_EQ((Q::coeff_at_t<0>::x::v), -1);
        EXPECT_EQ((Q::coeff_at_t<0>::y::v), 2);
        EXPECT_EQ((Q::coeff_at_t<1>::x::v), 1);
        EXPECT_EQ((Q::coeff_at_t<1>::y::v), 2);
    }
    // float divisibility
    {
        // x2 -1
        using A = polynomial<q32>::val<q32::one, q32::one>;
        // 2x + 2
        using B = polynomial<q32>::val<q32::val<i32::val<2>, i32::one>, q32::val<i32::val<2>, i32::one>>;
        using Q = polynomial<q32>::div_t<A, B>;
        EXPECT_EQ(Q::degree, 0);
        EXPECT_EQ((Q::coeff_at_t<0>::x::v), 1);
        EXPECT_EQ((Q::coeff_at_t<0>::y::v), 2);
    }
}

TEST(polynomials, gcd) {
    {
        // (x+1)*(x+1)
        using A = polynomial<q32>::val<q32::one, q32::val<i32::val<2>, i32::val<1>>, q32::one>;
        // (x+1)
        using B = polynomial<q32>::val<q32::one, q32::one>;
        using G = internal::gcd<polynomial<q32>>::type<A, B>;
        EXPECT_EQ(G::degree, 1);
        EXPECT_EQ((G::coeff_at_t<0>::x::v), 1);
        EXPECT_EQ((G::coeff_at_t<0>::y::v), 1);
        EXPECT_EQ((G::coeff_at_t<1>::x::v), 1);
        EXPECT_EQ((G::coeff_at_t<1>::y::v), 1);
    }
    {
        // (x+1)*(x+1)
        using A = polynomial<q32>::val<q32::one, q32::val<i32::val<2>, i32::val<1>>, q32::one>;
        // (x+1)*(x-1) :: x^2 - 1
        using B = polynomial<q32>::val<q32::one, q32::zero, q32::val<i32::val<-1>, i32::val<1>>>;
        // x + 1
        using G = polynomial<q32>::gcd_t<A, B>;
        EXPECT_EQ(G::degree, 1);
        EXPECT_EQ((G::coeff_at_t<0>::x::v), 1);
        EXPECT_EQ((G::coeff_at_t<0>::y::v), 1);
        EXPECT_EQ((G::coeff_at_t<1>::x::v), 1);
        EXPECT_EQ((G::coeff_at_t<1>::y::v), 1);
    }
}

TEST(polynomials, to_string) {
    using P32 = polynomial<q32>;
    using A = fpq32::val<P32::val<q32::one, q32::one>, P32::val<q32::val<i32::val<2>, i32::val<1>>, q32::one>>;
    auto rep = A::to_string();
    const char* expected = "(x + 1) / (2 x + 1)";
    EXPECT_STREQ(rep.c_str(), expected);
}

TEST(polynomials, order_over_rationals) {
    // constants compare by value, which polynomial gcds over q32 rely on
    using half = pq64::val<make_q64_t<1, 2>>;
    using third = pq64::val<make_q64_t<1, 3>>;
    EXPECT_TRUE((pq64::gt_t<half, third>::value));
    EXPECT_TRUE((pq64::gt_t<half, pq64::zero>::value));
    EXPECT_FALSE((pq64::gt_t<pq64::zero, half>::value));
}
