#include <gtest/gtest.h>

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

TEST(quotient_ring, instanciate_large) {
    using F17 = zpz<17>;
    using PF17 = polynomial<F17>;
    using F_17_8 = Quotient<PF17, ConwayPolynomial<17, 8>::type>;
    #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
    using x = F_17_8::inject_ring_t<PF17::X>;
    #pragma GCC diagnostic warning "-Wunused-local-typedefs"
}
