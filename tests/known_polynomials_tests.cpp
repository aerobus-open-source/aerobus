#include <gtest/gtest.h>
#include <type_traits>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

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

TEST(known_polynomials, bessel) {
    using B5 = known_polynomials::bessel<5>;
    EXPECT_EQ(5, B5::degree);
    EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<0>, i64::one>));
    EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<1>, i64::val<15>>));
    EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<2>, i64::val<105>>));
    EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<3>, i64::val<420>>));
    EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<4>, i64::val<945>>));
    EXPECT_TRUE((std::is_same_v<B5::template coeff_at_t<5>, i64::val<945>>));
}

TEST(known_polynomials, touchard) {
    using T5 = known_polynomials::touchard<5>;
    EXPECT_EQ(5, T5::degree);
    EXPECT_TRUE((std::is_same_v<T5::template coeff_at_t<0>, i64::zero>));
    EXPECT_TRUE((std::is_same_v<T5::template coeff_at_t<1>, i64::val<1>>));
    EXPECT_TRUE((std::is_same_v<T5::template coeff_at_t<2>, i64::val<15>>));
    EXPECT_TRUE((std::is_same_v<T5::template coeff_at_t<3>, i64::val<25>>));
    EXPECT_TRUE((std::is_same_v<T5::template coeff_at_t<4>, i64::val<10>>));
    EXPECT_TRUE((std::is_same_v<T5::template coeff_at_t<5>, i64::val<1>>));
    EXPECT_EQ((T5::eval<int64_t>(1)), (bell_v<i64, 5>));
}

TEST(known_polynomials, abel) {
    using A32 = known_polynomials::abel<3, 2>;
    EXPECT_EQ(3, A32::degree);
    EXPECT_TRUE((std::is_same_v<A32::template coeff_at_t<0>, i64::zero>));
    EXPECT_TRUE((std::is_same_v<A32::template coeff_at_t<1>, i64::val<36>>));
    EXPECT_TRUE((std::is_same_v<A32::template coeff_at_t<2>, i64::val<-12>>));
    EXPECT_TRUE((std::is_same_v<A32::template coeff_at_t<3>, i64::val<1>>));
}
