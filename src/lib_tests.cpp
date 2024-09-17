#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <gtest/gtest.h>
#include <cstdio>
#include <typeinfo>
#include <array>
#include <cmath>
#include <cstdlib>

#define AEROBUS_CONWAY_IMPORTS
#include "./aerobus.h"


// conformance : https://godbolt.org/z/Pn6s841bj

using namespace aerobus;  // NOLINT

TEST(utilies, type_at) {
    EXPECT_TRUE((std::is_same_v<internal::type_at_t<0, float, int, int64_t>, float>));
    EXPECT_TRUE((std::is_same_v<internal::type_at_t<1, float, int, int64_t>, int>));
    EXPECT_TRUE((std::is_same_v<internal::type_at_t<2, float, int, int64_t>, int64_t>));
}

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

TEST(polynomials, eval) {
    // 1 + 2x + 3x^2
    using poly = polynomial<i32>::val<i32::val<3>, i32::val<2>, i32::val<1>>;
    constexpr int v = poly::eval(1);
    EXPECT_EQ(v, 6);
    constexpr float vv = poly::eval(1.0f);
    EXPECT_EQ(vv, 6.0f);

    // 1/2 + 3x/2
    using polyf = polynomial<q32>::val<q32::val<i32::val<3>, i32::val<2>>, q32::val<i32::val<1>, i32::val<2>>>;
    constexpr float vvv = polyf::eval(1.0f);
    EXPECT_EQ(vvv, 2.0f);

    constexpr double vvvv = polyf::eval(-1.0);
    EXPECT_EQ(vvvv, -1.0);
}

TEST(fraction_field, eval) {
    using half = q32::val<i32::one, i32::val<2>>;
    constexpr float x = half::eval(2.0f);
    EXPECT_EQ(x, 0.5f);

    using thirdhalf = q32::val<i32::val<3>, i32::val<2>>;
    constexpr float y = thirdhalf::eval(1.0f);
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

TEST(i32, eq) {
    using a = i32::inject_constant_t<2>;
    using b = i32::inject_constant_t<3>;
    using c = i32::val<2>;
    EXPECT_TRUE((i32::eq_v<a, c>));
    EXPECT_FALSE((i32::eq_v<a, b>));
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

TEST(integers, gcd) {
    using A = i32::gcd_t<i32::val<12>, i32::val<6>>;
    EXPECT_EQ(A::v, 6);

    using B = i32::gcd_t<i32::val<15>, i32::val<6>>;
    EXPECT_EQ(B::v, 3);

    using C = i32::gcd_t<i32::val<5>, i32::val<3>>;
    EXPECT_EQ(C::v, 1);
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

TEST(quotient_ring, i32_z2z) {
    using QQ = Quotient<i32, i32::inject_constant_t<2>>;
    using two = QQ::inject_constant_t<2>;
    using three = QQ::inject_constant_t<3>;
    using one = QQ::inject_constant_t<1>;
    using zero = QQ::zero;
    EXPECT_TRUE((QQ::eq_v<two, zero>));
    EXPECT_TRUE((QQ::eq_v<one, three>));
}

TEST(quotient_ring, F4) {
    using F2 = zpz<2>;
    using PF2 = polynomial<F2>;
    using F4 = Quotient<PF2, ConwayPolynomial<2, 2>::type>;
    using zero = F4::zero;
    using one = F4::one;
    using phi = F4::inject_ring_t<PF2::X>;
    using phi2 = F4::inject_ring_t<PF2::template val<F2::inject_constant_t<1>, F2::inject_constant_t<1>>>;
    // check that elements are different
    EXPECT_FALSE((F4::eq_v<one, zero>));
    EXPECT_FALSE((F4::eq_v<phi, zero>));
    EXPECT_FALSE((F4::eq_v<phi2, zero>));
    EXPECT_FALSE((F4::eq_v<phi, one>));
    EXPECT_FALSE((F4::eq_v<phi2, one>));
    EXPECT_FALSE((F4::eq_v<phi, phi2>));

    // 1 + 1 = 0
    EXPECT_TRUE((F4::eq_v<zero, F4::template add_t<one, one>>));
    // one + phi = phi2
    EXPECT_TRUE((F4::eq_v<phi2, F4::template add_t<one, phi>>));
    // one + phi2 = phi
    EXPECT_TRUE((F4::eq_v<phi, F4::template add_t<one, phi2>>));
    // phi * phi = phi2
    EXPECT_TRUE((F4::eq_v<phi2, F4::template mul_t<phi, phi>>));
    // phi2 * phi2 = phi
    EXPECT_TRUE((F4::eq_v<phi, F4::template mul_t<phi2, phi2>>));
}

TEST(quotient_ring, instanciate_large) {
    using F17 = zpz<17>;
    using PF17 = polynomial<F17>;
    using F_17_8 = Quotient<PF17, ConwayPolynomial<17, 8>::type>;
    using x = F_17_8::inject_ring_t<PF17::X>;
}

TEST(utilities, factorial) {
    constexpr float x = factorial_v<i32, 3>;
    EXPECT_EQ(x, 6.0f);

    constexpr size_t y = factorial_v<i32, 0>;
    EXPECT_EQ(y, 1);
}

TEST(utilities, combination) {
    constexpr int x = combination_v<i32, 2, 4>;
    EXPECT_EQ(x, 6);

    constexpr int y = combination_v<i32, 0, 4>;
    EXPECT_EQ(y, 1);

    constexpr int z = combination_v<i32, 1, 4>;
    EXPECT_EQ(z, 4);

    constexpr int zz = combination_v<i32, 3, 4>;
    EXPECT_EQ(zz, 4);
}

TEST(utilities, bernoulli) {
    constexpr float b0 = bernoulli_v<float, i32, 0>;
    EXPECT_EQ(b0, 1.0f);

    constexpr float b1 = bernoulli_v<float, i32, 1>;
    EXPECT_EQ(b1, -0.5f);

    using B2 = bernoulli_t<i32, 2>;
    EXPECT_EQ(B2::x::v, 1);
    EXPECT_EQ(B2::y::v, 6);

    constexpr double b3 = bernoulli_v<double, i32, 3>;
    EXPECT_EQ(b3, 0.0);

    using B4 = bernoulli_t<i64, 4>;
    EXPECT_EQ(B4::x::v, -1);
    EXPECT_EQ(B4::y::v, 30);
}

TEST(utilities, bell) {
    constexpr int64_t B0 = bell_v<i64, 0>;
    EXPECT_EQ(B0, 1);
    constexpr int64_t B1 = bell_v<i64, 1>;
    EXPECT_EQ(B1, 1);
    constexpr int64_t B2 = bell_v<i64, 2>;
    EXPECT_EQ(B2, 2);
    constexpr int64_t B3 = bell_v<i64, 3>;
    EXPECT_EQ(B3, 5);
    constexpr int64_t B4 = bell_v<i64, 4>;
    EXPECT_EQ(B4, 15);
}

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

TEST(utilities, is_prime) {
    EXPECT_TRUE(is_prime_v<2>);
    EXPECT_TRUE(is_prime_v<3>);
    EXPECT_TRUE(is_prime_v<7>);
    EXPECT_TRUE(is_prime_v<11>);
    EXPECT_TRUE(is_prime_v<13>);
    EXPECT_TRUE(is_prime_v<31>);
    EXPECT_TRUE(is_prime_v<7883>);
    EXPECT_TRUE(is_prime_v<7919>);
    EXPECT_TRUE(is_prime_v<7927>);
    EXPECT_TRUE(is_prime_v<1000423>);
    EXPECT_FALSE(is_prime_v<1>);
    EXPECT_FALSE(is_prime_v<4>);
    EXPECT_FALSE(is_prime_v<7884>);
    EXPECT_FALSE(is_prime_v<7928>);
}

TEST(utilities, exp) {
    using E = aerobus::exp<i32, 12>;
    constexpr float e0 = E::eval(0.0F);
    constexpr float e01 = E::eval(0.1F);
    EXPECT_EQ(e0, 1.0f);

    EXPECT_TRUE((std::abs(std::exp(0.1f) - e01) <= 1E-7F));
}

TEST(utilities, alternate) {
    constexpr int a0 = internal::alternate<i32, 0>::value;
    EXPECT_EQ(a0, 1);
    constexpr int a1 = internal::alternate<i32, 1>::value;
    EXPECT_EQ(a1, -1);
    constexpr int a2 = internal::alternate<i32, 2>::value;
    EXPECT_EQ(a2, 1);
}

TEST(utilities, pow) {
    constexpr int32_t a0 = pow_t<i32, typename i32::inject_constant_t<2>, 3>::v;
    EXPECT_EQ(a0, 8);
    constexpr int32_t a1 = pow_t<i32, typename i32::inject_constant_t<2>, 4>::v;
    EXPECT_EQ(a1, 16);
    constexpr int32_t a2 = pow_t<i32, i32::zero, 0>::v;
    EXPECT_EQ(a2, 1);
    constexpr int64_t a3 = pow_v<i64, typename i64::inject_constant_t<3>, 3>;
    EXPECT_EQ(a3, 27);
}

TEST(utilities, pow_scalar) {
    constexpr int32_t a0 = pow_scalar<int32_t, 2>(2);
    EXPECT_EQ(a0, 4);

    constexpr double x = pow_scalar<double, 3>(3.0);
    EXPECT_EQ(x, 27.0);

    constexpr double one = pow_scalar<size_t, 0>(12);
    EXPECT_EQ(one, 1);

    constexpr double zone = pow_scalar<size_t, 0>(0);
    EXPECT_EQ(zone, 1);
}

TEST(utilities, abs) {
    using a0 = abs_t<i32::inject_constant_t<-1>>;
    EXPECT_EQ(a0::v, 1);
    using a1 = abs_t<i64::inject_constant_t<1>>;
    EXPECT_EQ(a1::v, 1);
    using a2 = abs_t<pi64::val<i64::val<-1>, i64::val<1>>>;
    using expected = pi64::val<i64::val<1>, i64::val<-1>>;
    EXPECT_TRUE((std::is_same_v<a2, expected>));
}

TEST(utilities, stirling) {
    constexpr int s00 = stirling_signed_v<i64, 0, 0>;
    EXPECT_EQ(s00, 1);
    constexpr int s01 = stirling_signed_v<i64, 0, 1>;
    EXPECT_EQ(s01, 0);
    constexpr int s11 = stirling_signed_v<i64, 1, 1>;
    EXPECT_EQ(s11, 1);
    constexpr int s53 = stirling_signed_v<i64, 5, 3>;
    EXPECT_EQ(s53, 35);
    constexpr int s94 = stirling_signed_v<i64, 9, 4>;
    EXPECT_EQ(s94, -67284);
    constexpr int u94 = stirling_unsigned_v<i64, 9, 4>;
    EXPECT_EQ(u94, 67284);
}

TEST(type_list, basic_assertions) {
    using A = type_list<int, float>;
    EXPECT_EQ(A::length, 2);
    EXPECT_TRUE((std::is_same_v<A::at<0>, int>));
    EXPECT_TRUE((std::is_same_v<A::at<1>, float>));

    using B = A::push_back<double>;
    EXPECT_EQ(B::length, 3);
    EXPECT_TRUE((std::is_same_v<B::at<2>, double>));

    using C = B::pop_front;
    EXPECT_EQ(C::tail::length, 2);
    EXPECT_TRUE((std::is_same_v<C::type, int>));
    EXPECT_TRUE((std::is_same_v<C::tail::at<0>, float>));

    using D = C::tail::split<1>;
    EXPECT_EQ(D::head::length, 1);
    EXPECT_EQ(D::tail::length, 1);
    EXPECT_TRUE((std::is_same_v<D::head::at<0>, float>));
    EXPECT_TRUE((std::is_same_v<D::tail::at<0>, double>));
}

TEST(type_list, more_complex_assertions) {
    {
        using t = type_list<float, int, double>;
        using t2 = t::remove<1>;
        EXPECT_EQ(t2::length, 2);
        EXPECT_TRUE((std::is_same_v<t2::at<0>, float>));
        EXPECT_TRUE((std::is_same_v<t2::at<1>, double>));
    }
    {
        using t = type_list<float, double>;
        using t2 = t::insert<int, 1>;
        EXPECT_EQ(t2::length, 3);
        EXPECT_TRUE((std::is_same_v<t2::at<0>, float>));
        EXPECT_TRUE((std::is_same_v<t2::at<1>, int>));
        EXPECT_TRUE((std::is_same_v<t2::at<2>, double>));
    }
    {
        using t = type_list<>;
        using t2 = t::insert<float, 0>;
        EXPECT_EQ(t2::length, 1);
        EXPECT_TRUE((std::is_same_v<t2::at<0>, float>));
    }
}

TEST(concepts, ring) {
    EXPECT_TRUE((aerobus::IsRing<i32>));
    EXPECT_TRUE((aerobus::IsRing<i64>));
    EXPECT_TRUE((aerobus::IsRing<zpz<3>>));
    EXPECT_TRUE((aerobus::IsRing<aerobus::q32>));
    EXPECT_TRUE((aerobus::IsField<aerobus::q32>));
    EXPECT_TRUE((aerobus::IsRing<aerobus::polynomial<i32>>));
    EXPECT_TRUE((aerobus::IsField<aerobus::FractionField<aerobus::polynomial<i32>>>));
}

TEST(continuous_fractions, basic_assertions) {
    // A001203
    constexpr double A_PI = PI_fraction::val;
    EXPECT_FALSE((::fabs(A_PI - M_PI) > 0));

    // A003417
    constexpr double A_E = E_fraction::val;
    EXPECT_FALSE((::fabs(A_E - M_E) > 0));

    constexpr double A_SQRT2 = SQRT2_fraction::val;
    EXPECT_FALSE((::fabs(A_SQRT2 - M_SQRT2) > 0));

    constexpr double A_SQRT3 = SQRT3_fraction::val;
    EXPECT_FALSE((::fabs(A_SQRT3 - 1.7320508075688772935) > 0));
}

TEST(known_polynomials, chebyshev) {
    using T4 = known_polynomials::chebyshev_T<4>;

    EXPECT_EQ(T4::degree, 4);
    EXPECT_EQ((T4::template coeff_at_t<4>::v), 8);
    EXPECT_EQ((T4::template coeff_at_t<3>::v), 0);
    EXPECT_EQ((T4::template coeff_at_t<2>::v), -8);
    EXPECT_EQ((T4::template coeff_at_t<1>::v), 0);
    EXPECT_EQ((T4::template coeff_at_t<0>::v), 1);

    using U4 = known_polynomials::chebyshev_U<4>;
    EXPECT_EQ(U4::degree, 4);
    EXPECT_EQ((U4::template coeff_at_t<4>::v), 16);
    EXPECT_EQ((U4::template coeff_at_t<3>::v), 0);
    EXPECT_EQ((U4::template coeff_at_t<2>::v), -12);
    EXPECT_EQ((U4::template coeff_at_t<1>::v), 0);
    EXPECT_EQ((U4::template coeff_at_t<0>::v), 1);
}


TEST(known_polynomials, laguerre) {
    {
        using L2 = known_polynomials::laguerre<2>;
        EXPECT_EQ((L2::coeff_at_t<2>::template get<float>()), 0.5F);
        EXPECT_EQ((L2::coeff_at_t<1>::template get<float>()),  -2.0F);
        EXPECT_EQ((L2::coeff_at_t<0>::template get<float>()),  1.0F);
    }
    {
        using L3 = known_polynomials::laguerre<3>;
        EXPECT_EQ((L3::coeff_at_t<3>::x::v), -1);
        EXPECT_EQ((L3::coeff_at_t<3>::y::v), 6);
        EXPECT_EQ((L3::coeff_at_t<2>::template get<float>()), 1.5F);
        EXPECT_EQ((L3::coeff_at_t<1>::template get<float>()), -3.0F);
        EXPECT_EQ((L3::coeff_at_t<0>::template get<float>()), 1.0F);
    }
}


TEST(known_polynomials, hermite) {
    {
        using H2 = known_polynomials::hermite_prob<2>;
        EXPECT_EQ((H2::coeff_at_t<2>::get<float>()), 1.0F);
        EXPECT_EQ((H2::coeff_at_t<1>::get<float>()), 0.0F);
        EXPECT_EQ((H2::coeff_at_t<0>::get<float>()), -1.0F);
    }
    {
        using H3 = known_polynomials::hermite_prob<3>;
        EXPECT_EQ((H3::coeff_at_t<3>::get<float>()), 1.0F);
        EXPECT_EQ((H3::coeff_at_t<2>::get<float>()), 0.0F);
        EXPECT_EQ((H3::coeff_at_t<1>::get<float>()), -3.0F);
        EXPECT_EQ((H3::coeff_at_t<0>::get<float>()), 0.0F);
    }
    {
        using H2 = known_polynomials::hermite_phys<2>;
        EXPECT_EQ((H2::coeff_at_t<2>::get<float>()), 4.0F);
        EXPECT_EQ((H2::coeff_at_t<1>::get<float>()), 0.0F);
        EXPECT_EQ((H2::coeff_at_t<0>::get<float>()), -2.0F);
    }
    {
        using H3 = known_polynomials::hermite_phys<3>;
        EXPECT_EQ((H3::coeff_at_t<3>::get<float>()), 8.0F);
        EXPECT_EQ((H3::coeff_at_t<2>::get<float>()), 0.0F);
        EXPECT_EQ((H3::coeff_at_t<1>::get<float>()), -12.0F);
        EXPECT_EQ((H3::coeff_at_t<0>::get<float>()), 0.0F);
    }
}

TEST(known_polynomials, bernstein) {
    {
        using B00 = known_polynomials::bernstein<0, 0>;
        EXPECT_EQ(B00::degree, 0);
        EXPECT_EQ(B00::coeff_at_t<0>::v, 1);
    }
    {
        // 1 - X
        using B01 = known_polynomials::bernstein<0, 1>;
        EXPECT_EQ(B01::degree, 1);
        EXPECT_EQ(B01::coeff_at_t<0>::v, 1) << "B01(0) != 1";
        EXPECT_EQ(B01::coeff_at_t<1>::v, -1) << "B01(1) != -1";
    }
    {
        // X
        using B11 = known_polynomials::bernstein<1, 1>;
        EXPECT_EQ(B11::degree, 1);
        EXPECT_EQ(B11::coeff_at_t<0>::v, 0);
        EXPECT_EQ(B11::coeff_at_t<1>::v, 1);
    }
    {
        // 2X(1-X)
        using B12 = known_polynomials::bernstein<1, 2>;
        EXPECT_EQ(B12::degree, 2);
        EXPECT_EQ(B12::coeff_at_t<0>::v, 0);
        EXPECT_EQ(B12::coeff_at_t<1>::v, 2);
        EXPECT_EQ(B12::coeff_at_t<2>::v, -2);
    }
    {
        // partition
        using B03 = known_polynomials::bernstein<0, 3>;
        using B13 = known_polynomials::bernstein<1, 3>;
        using B23 = known_polynomials::bernstein<2, 3>;
        using B33 = known_polynomials::bernstein<3, 3>;
        using sum = vadd_t<B03, B13, B23, B33>;
        EXPECT_TRUE((std::is_same_v<sum, typename pi64::one>));
    }
}

TEST(known_polynomials, legendre) {
    {
        using L1 = known_polynomials::legendre<1>;
        EXPECT_TRUE((std::is_same_v<L1, typename pq64::X>));
    }
    {
        using L5 = known_polynomials::legendre<5>;
        EXPECT_TRUE((std::is_same_v<
            L5::template coeff_at_t<3>,
            make_q64_t<-70, 8>>)) << L5::template coeff_at_t<3>::to_string();
    }
    {
        using L10 = known_polynomials::legendre<10>;
        EXPECT_TRUE((std::is_same_v<
            L10::template coeff_at_t<8>,
            make_q64_t<-109395, 256>>)) << L10::template coeff_at_t<8>::to_string();
    }
}

TEST(known_polynomials, bernoulli) {
    {
        using B1 = known_polynomials::bernoulli<1>;
        EXPECT_EQ(1, B1::degree);
        EXPECT_TRUE((std::is_same_v<B1::template coeff_at_t<0>, make_q64_t<-1, 2>>));
        EXPECT_TRUE((std::is_same_v<B1::template coeff_at_t<1>, typename q64::one>));
    } {
        using B5 = known_polynomials::bernoulli<5>;
        EXPECT_EQ(5, B5::degree);
        EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<0>, q64::zero>));
        EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<1>, make_q64_t<-1, 6>>));
        EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<2>, q64::zero>));
        EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<3>, make_q64_t<5, 3>>));
        EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<4>, make_q64_t<-5, 2>>));
        EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<5>, q64::one>));
    }
}

TEST(embedding, simple) {
    {
        using Small = Quotient<i32, i32::val<4>>;
        using embedded = typename Embed<Small, i32>::type<typename Small::template inject_constant_t<3>>;
        EXPECT_TRUE((std::is_same_v<embedded, typename i32::template inject_constant_t<3>>));
    }
    {
        using Small = zpz<4>;
        using embedded = typename Embed<Small, i32>::type<typename Small::template inject_constant_t<3>>;
        EXPECT_EQ(embedded::v, 3);
    }
    {
        using embedded = typename Embed<i32, i64>::type<i32::val<12>>;
        EXPECT_TRUE((std::is_same_v<embedded, typename i64::template inject_constant_t<12>>));
    }
    {
        using embedded = typename Embed<i32, q32>::type<i32::val<12>>;
        EXPECT_TRUE((std::is_same_v<embedded, q32::inject_constant_t<12>>));
    }
    {
        using embedded = typename Embed<q32, q64>::type<make_q32_t<1, 2>>;
        EXPECT_TRUE((std::is_same_v<embedded, make_q64_t<1, 2>>));
    }
}

TEST(embedding, polynomial) {
    {  // 1+X
        using embedded = Embed<polynomial<i32>, polynomial<i64>>::type<polynomial<i32>::val<i32::one, i32::one>>;
        EXPECT_TRUE((std::is_same_v<embedded, polynomial<i64>::val<i64::one, i64::one>>));
    }
    {  // 1+2X
        using embedded = Embed<polynomial<i32>, polynomial<i64>>::type<polynomial<i32>::val<
                            i32::inject_constant_t<2>, i32::one>>;
        EXPECT_TRUE((std::is_same_v<embedded, polynomial<i64>::val<i64::inject_constant_t<2>, i64::one>>));
    }
    {  // 1+2X, from i32 to q32
        using embedded = embed_int_poly_in_fractions_t<make_int_polynomial_t<i32, 2, 1>>;
        EXPECT_TRUE((std::is_same_v<embedded, make_frac_polynomial_t<i32, 2, 1>>))
            << "actual : " << embedded::to_string();
    }
}
