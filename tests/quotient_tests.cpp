#include <gtest/gtest.h>
#include <type_traits>

#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"

using namespace aerobus;  // NOLINT

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

TEST(quotient_ring, F9) {
    using F3 = zpz<3>;
    using PF3 = polynomial<F3>;
    using F9 = Quotient<PF3, ConwayPolynomial<3, 2>::type>;  // x^2 + 2x + 2
    using phi = F9::inject_ring_t<PF3::X>;
    using phi2 = F9::mul_t<phi, phi>;
    using phi4 = F9::mul_t<phi2, phi2>;
    using phi8 = F9::mul_t<phi4, phi4>;
    // phi^2 = -2 phi - 2 = phi + 1
    EXPECT_TRUE((F9::eq_v<phi2, F9::inject_ring_t<PF3::val<F3::one, F3::one>>>));
    // a Conway polynomial is primitive: phi has order 8
    EXPECT_TRUE((F9::eq_v<phi4, F9::inject_constant_t<-1>>));
    EXPECT_TRUE((F9::eq_v<phi8, F9::one>));
}

TEST(quotient_ring, F25) {
    using F5 = zpz<5>;
    using PF5 = polynomial<F5>;
    using F25 = Quotient<PF5, ConwayPolynomial<5, 2>::type>;  // x^2 + 4x + 2
    using phi = F25::inject_ring_t<PF5::X>;
    using phi2 = F25::mul_t<phi, phi>;
    using phi4 = F25::mul_t<phi2, phi2>;
    using phi8 = F25::mul_t<phi4, phi4>;
    using phi12 = F25::mul_t<phi8, phi4>;
    using phi24 = F25::mul_t<phi12, phi12>;
    // phi^2 = -4 phi - 2 = phi + 3
    EXPECT_TRUE((F25::eq_v<phi2, F25::inject_ring_t<PF5::val<F5::one, F5::val<3>>>>));
    // phi has order 24
    EXPECT_FALSE((F25::eq_v<phi8, F25::one>));
    EXPECT_TRUE((F25::eq_v<phi12, F25::inject_constant_t<-1>>));
    EXPECT_TRUE((F25::eq_v<phi24, F25::one>));
}

TEST(quotient_ring, instanciate_large) {
    using F17 = zpz<17>;
    using PF17 = polynomial<F17>;
    using F_17_8 = Quotient<PF17, ConwayPolynomial<17, 8>::type>;
    #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
    using x = F_17_8::inject_ring_t<PF17::X>;
    #pragma GCC diagnostic warning "-Wunused-local-typedefs"
}

TEST(quotient_ring, integers_negative_representatives) {
    using Z17 = Quotient<i64, i64::val<17>>;
    // the class of -1 is the class of 16
    EXPECT_TRUE((Z17::eq_v<Z17::inject_constant_t<-1>, Z17::inject_constant_t<16>>));
    EXPECT_TRUE((std::is_same_v<Z17::inject_constant_t<-1>, Z17::inject_constant_t<16>>));
    EXPECT_TRUE((Z17::eq_v<Z17::add_t<Z17::inject_constant_t<-1>, Z17::one>, Z17::zero>));
    EXPECT_EQ((Z17::inject_constant_t<-1>::type::v), 16);
}

TEST(quotient_ring, integer_quotient_is_a_ring) {
    using Z17 = Quotient<i64, i64::val<17>>;
    using two = Z17::inject_constant_t<2>;
    using nine = Z17::inject_constant_t<9>;
    // 2 * 9 = 18 = 1
    EXPECT_TRUE((std::is_same_v<Z17::div_t<Z17::one, two>, nine>));
    EXPECT_TRUE((std::is_same_v<Z17::sub_t<Z17::zero, Z17::one>, Z17::inject_constant_t<16>>));
    EXPECT_TRUE((std::is_same_v<Z17::mod_t<nine, two>, Z17::zero>));
    EXPECT_TRUE((std::is_same_v<Z17::gcd_t<two, nine>, Z17::one>));
    EXPECT_TRUE((Z17::gt_v<nine, two>));
    EXPECT_TRUE((Z17::lt_v<two, nine>));
    EXPECT_TRUE(Z17::is_field);
    EXPECT_FALSE((Quotient<i64, i64::val<16>>::is_field));
    EXPECT_TRUE((IsEuclideanDomain<Z17>));
    EXPECT_TRUE((IsField<Z17>));
    // values know their ring
    EXPECT_TRUE((std::is_same_v<aerobus::add_t<two, nine>, Z17::inject_constant_t<11>>));
    EXPECT_EQ(two::to_string(), "2");
    EXPECT_TRUE(Z17::zero::is_zero_v);
    EXPECT_EQ((nine::get<int>()), 9);
}

TEST(quotient_ring, field_division_F4) {
    using F2 = zpz<2>;
    using PF2 = polynomial<F2>;
    using F4 = Quotient<PF2, ConwayPolynomial<2, 2>::type, true>;
    using phi = F4::inject_ring_t<PF2::X>;
    using phi2 = F4::mul_t<phi, phi>;
    EXPECT_TRUE(F4::is_field);
    EXPECT_TRUE((IsField<F4>));
    // phi^3 = 1, so 1 / phi = phi^2 and phi^2 / phi = phi
    EXPECT_TRUE((std::is_same_v<F4::div_t<F4::one, phi>, phi2>));
    EXPECT_TRUE((std::is_same_v<F4::div_t<phi2, phi>, phi>));
    EXPECT_TRUE((std::is_same_v<F4::mod_t<phi2, phi>, F4::zero>));
}

TEST(quotient_ring, polynomials_over_F4) {
    using F2 = zpz<2>;
    using PF2 = polynomial<F2>;
    using F4 = Quotient<PF2, ConwayPolynomial<2, 2>::type, true>;
    using phi = F4::inject_ring_t<PF2::X>;
    using phi2 = F4::mul_t<phi, phi>;
    using PF4 = polynomial<F4>;
    // (y + phi)(y + phi^2) = y^2 + (phi + phi^2) y + phi^3 = y^2 + y + 1
    using a = PF4::val<F4::one, phi>;
    using b = PF4::val<F4::one, phi2>;
    using product = PF4::val<F4::one, F4::one, F4::one>;
    EXPECT_TRUE((std::is_same_v<PF4::mul_t<a, b>, product>));
    EXPECT_TRUE((std::is_same_v<PF4::div_t<product, a>, b>));
}

TEST(quotient_ring, galois_field) {
    using F9 = GaloisField<3, 2>;
    EXPECT_TRUE(F9::is_field);
    using phi = F9::inject_ring_t<polynomial<zpz<3>>::X>;
    // phi * (1 / phi) = 1
    EXPECT_TRUE((std::is_same_v<F9::mul_t<phi, F9::div_t<F9::one, phi>>, F9::one>));
}

TEST(quotient_ring, bigint_integers) {
    // remainders of bigint have the sign of the dividend, as for i64: classes need canonical representatives
    using Z17 = Quotient<bigint, bigint::inject_constant_t<17>>;
    static_assert(std::is_same_v<Z17::inject_constant_t<-1>, Z17::inject_constant_t<16>>);
    static_assert(Z17::is_field);
    static_assert(IsFiniteField<Z17> && field_order_v<Z17> == 17);
    static_assert(std::is_same_v<Z17::div_t<Z17::one, Z17::inject_constant_t<2>>, Z17::inject_constant_t<9>>);
    // a negative, composite modulus
    using Z15 = Quotient<bigint, bigint::inject_constant_t<-15>>;
    static_assert(!Z15::is_field);
    static_assert(std::is_same_v<Z15::inject_constant_t<-1>, Z15::inject_constant_t<14>>);
    // beyond 64 bits primality is not tested: 2^127 - 1 is declared prime
    using M127 = bigint::sub_t<pow_t<bigint::inject_constant_t<2>, 127>, bigint::one>;
    using Q127 = Quotient<bigint, M127, true>;
    static_assert(Q127::is_field);
    static_assert(std::is_same_v<Q127::inject_constant_t<-1>, Q127::val<bigint::sub_t<M127, bigint::one>>>);
}

TEST(quotient_ring, large_prime_moduli) {
    // 2^61 - 1 is prime, too large for trial division at compile time; 2^61 + 1 is a multiple of 3
    static_assert(Quotient<i64, i64::val<2305843009213693951>>::is_field);
    static_assert(!Quotient<i64, i64::val<2305843009213693953>>::is_field);
    static_assert(Quotient<bigint, bigint::from_string_t<"2305843009213693951">>::is_field);
}
