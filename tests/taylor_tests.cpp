#include <gtest/gtest.h>
#include <cmath>

#include "aerobus.h"

using namespace aerobus;  // NOLINT

TEST(utilities, exp) {
    using E = aerobus::exp<i32, 12>;
    constexpr float e0 = E::eval(0.0F);
    constexpr float e01 = E::eval(0.1F);
    EXPECT_EQ(e0, 1.0f);

    EXPECT_TRUE((std::abs(std::exp(0.1f) - e01) <= 1E-7F));
}

TEST(utilities, expm1) {
    using E = aerobus::expm1<i32, 12>;
    constexpr float e0 = E::eval(0.0F);
    constexpr float e01 = E::eval(0.1F);
    EXPECT_EQ(e0, 0.0f);

    EXPECT_TRUE((std::abs(std::expm1(0.1f) - e01) <= 1E-7F));
}

TEST(utilities, expm1_long_double) {
    using E = aerobus::expm1<i64, 12>;
    long double e0 = E::compensated_eval(0.0L);
    long double e01 = E::compensated_eval(0.1L);
    EXPECT_EQ(e0, 0.0);

    EXPECT_TRUE((std::abs(std::expm1(0.1) - e01) <= 1E-7));
}

TEST(sollya, expm1_float) {
    using P = polynomial<q32>::template val<
        make_q32_t<33771, 1320736876>,
        make_q32_t<192347, 939654189>,
        make_q32_t<2498407, 1799830379>,
        make_q32_t<7975888, 957650921>,
        make_q32_t<34034752, 816828017>,
        make_q32_t<344868026, 2069191783>,
        make_q32_t<326970344, 653940747>,
        make_q32_t<2009304618, 2009304817>,
        make_q32_t<1, 907978194>>;
    float infnorm = 0.0F;
    for (float x = -0.099F; x <= 0.099F; x += 0.01F) {
        float y = P::eval(x);
        float xy = std::abs(y - std::expm1(x));
        if (xy > infnorm) {
            infnorm = xy;
        }
    }

    EXPECT_TRUE((infnorm < 1.5E-8F)) << infnorm;
}

TEST(sollya, taylor_vs_remez) {
    using P = aerobus::polynomial<aerobus::q32>::template val<
        aerobus::make_q32_t<33771, 1320736876>,
        aerobus::make_q32_t<192347, 939654189>,
        aerobus::make_q32_t<2498407, 1799830379>,
        aerobus::make_q32_t<7975888, 957650921>,
        aerobus::make_q32_t<34034752, 816828017>,
        aerobus::make_q32_t<344868026, 2069191783>,
        aerobus::make_q32_t<326970344, 653940747>,
        aerobus::make_q32_t<2009304618, 2009304817>,
        aerobus::make_q32_t<1, 907978194>>;

    float inf_remez = 0.0F;
    float inf_taylor = 0.0F;

    for (float x = -0.99F; x <= 0.99F; x += 0.01F) {
        float y = P::eval(x);
        float yy = aerobus::expm1<i32, 8>::eval(x);
        float remez_error = std::abs(y - std::expm1(x));
        float taylor_error = std::abs(yy - std::expm1(x));
        if (remez_error > inf_remez) {
            inf_remez = remez_error;
        }
        if (taylor_error > inf_taylor) {
            inf_taylor = taylor_error;
        }
    }

    EXPECT_TRUE((inf_remez < inf_taylor));
}
