#include <iostream>
#include <cstdio>
#include <chrono>
#include "./aerobus.h"


template<typename FT>
std::enable_if_t<std::is_floating_point_v<FT>, FT>
ulp(FT x) {
    if (x > 0) {
        return std::nexttoward(x, std::numeric_limits<FT>::infinity()) - x;
    } else {
        return x - std::nexttoward(x, -std::numeric_limits<FT>::infinity());
    }
}


// TEST(libm, exp2) {
//     EXPECT_EQ(aerobus::libm::exp2(1.0), 2.0);
//     EXPECT_EQ(aerobus::libm::exp2(1.00000000001), std::exp2(1.00000000001));
//     EXPECT_EQ(aerobus::libm::exp2(2.0), 4.0);
//     EXPECT_EQ(aerobus::libm::exp2(2.0), 4.0);
// }

TEST(libm, sin_reduction) {
    using d_constants = aerobus::arithmetic_helpers<double>;
    using f_constants = aerobus::arithmetic_helpers<float>;
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(0.1, 0.1F);
        EXPECT_TRUE(red.return_fast_sin);
        EXPECT_EQ(red.transform, 0.1) << "input was : " << 0.1 << std::endl;
    }
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(1.0, 1.0F);
        EXPECT_TRUE(red.return_fast_cos);
        EXPECT_EQ(red.transform, d_constants::pi_2() - 1.0) << "input was : " << 1.0 << std::endl;
        EXPECT_FALSE(red.negate);
    }
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(2.35, 2.35F);
        EXPECT_TRUE(red.return_fast_cos);
        EXPECT_EQ(
                red.transform,
                d_constants::pi_2() - (d_constants::pi() - 2.35))
            << "input was : " << 2.35 << std::endl;
    }
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(3.15, 3.15F);
        EXPECT_TRUE(red.return_fast_sin);
        EXPECT_EQ(
                red.transform,
                d_constants::pi() - (d_constants::two_pi() - 3.15))
            << "input was : " << 2.35 << std::endl;
        EXPECT_TRUE(red.negate);
    }
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(d_constants::pi(), f_constants::pi());
        EXPECT_TRUE(red.return_x);
        EXPECT_EQ(red.transform, 0) << "input was : " << d_constants::pi() << std::endl;
        EXPECT_TRUE(red.negate);
    }
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(-1.0, -1.0F);
        EXPECT_TRUE(red.return_fast_cos);
        EXPECT_EQ(red.transform, d_constants::pi_2() - 1.0) << "input was : " << -1.0 << std::endl;
        EXPECT_TRUE(red.negate);
    }
    {
        auto red = aerobus::libm::sin_reduction<float>::eval(-0.0, -0.0F);
        EXPECT_TRUE(red.return_x);
        EXPECT_EQ(red.transform, static_cast<double>(-0.0F)) << "input was : " << -0.0 << std::endl;
        EXPECT_FALSE(red.negate);
    }
}

template<typename FT>
std::enable_if_t<std::is_floating_point_v<FT>, FT>
ulp(FT x) {
    if (x > 0) {
        return std::nexttoward(x, std::numeric_limits<FT>::infinity()) - x;
    } else {
        return x - std::nexttoward(x, -std::numeric_limits<FT>::infinity());
    }
}

time_t time_since_epoch() {
    auto now = std::chrono::system_clock::now();
    return std::chrono::system_clock::to_time_t( now );
}

int test_sin_exhaustive() {
    time_t current = time_since_epoch();
    float start = 2.818502783775e-01;  // std::numeric_limits<float>::min()
    uint64_t okcount = 0, failcount = 0;
    uint64_t count = 0;
    float err = 0;
    bool firstlog = true;
    for (float x = start; x <= std::numeric_limits<float>::max();) {
        float s = aerobus::libm::sin(x);
        float expected = std::sin(x);
        time_t t = time_since_epoch();
        if (t != current) {
            printf("%.12e, ok: %ld - fail: %ld, percentage error: %.12f\n", x, okcount, failcount, 100.0F * (float) failcount / (float) (failcount + okcount));
            current = t;
        }
        if (std::isinf(x)) {  // should not happen given doc of std::numeric_limits but still double check
            if (s == s) {
                return 1;
            }
        } else {
            if (s < expected - ulp(expected) || s > expected + ulp(expected)) {
                failcount += 1;
                err += std::fabs(s - expected);
                if(firstlog) {
                    std::cout << std::hexfloat << "input : " << x << " expected : " << expected << " computed " << s <<
                        std::endl << "    difference is : " << s - expected << std::endl;
                    firstlog = false;
                }
                // return 1;
            } else {
                okcount += 1;
            }
        }

        x = std::nexttoward(x, std::numeric_limits<float>::infinity());
    }

    printf("ok : %ld, fail: %ld\n", okcount, failcount);

    return 0;
}


TEST(libm, sin) {
    using constants = aerobus::arithmetic_helpers<float>;
    float values[] = {
        0.0000001F,
        -0.00000001F,
        0.863769531250000F,
        0.861816406250000F,
        constants::pi() / 8,
        3 * constants::pi() / 8,
        3 * constants::pi() / 4,
        3 * constants::pi() / 2,
        constants::pi(),
        constants::pi_2(),
        constants::two_pi(),
        9 * constants::pi() / 8,
        10.0F,
        -1.0F,
        -6.050781250000F,
        776.0F,
    };

    for (float x : values) {
        float aero = aerobus::libm::sin(x);
        float expected = static_cast<float>(std::sin(static_cast<double>(x)));
        EXPECT_TRUE(aero <= expected + ulp(expected) && aero >= expected - ulp(expected)) <<
            std::hexfloat << "input : " << x << " expected : " << expected << " computed " << aero <<
            std::endl << "    difference is : " << aero - expected << std::endl;
        // uncomment to see small differences
        // EXPECT_EQ(aero, expected) <<
        //     std::hexfloat << "input : " << x << " expected : " << expected << " computed " << aero <<
        //     std::endl << "    difference is : " << aero - expected << std::endl;
    }

    float exact_values[14] = {
        NAN, NAN,
        1E-26F, 1E-26F,
        -1E-26F, -1E-26F,
        -0.F, -0.F,
        0.F, 0.F,
        // 0x1.4f1a6ep+1, 0x1.fffff4p-2,  // given by sollya and clang libm where the agree
        std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN(),
        -std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN()
    };

    for (unsigned long i = 0; i < sizeof(exact_values) / sizeof(float); i += 2) {
        float x = exact_values[i];
        float aero = aerobus::libm::sin(x);
        float expected = exact_values[i+1];
        if (x != x || expected != expected) {
            EXPECT_TRUE(aero != aero);
        } else {
            EXPECT_EQ(aero, expected) << "input was : " << x
                << std::hexfloat << " (" << x << ") and yielded : "
                << aero << " instead of " << expected << std::endl;
        }
    }
}

TEST(libm, cos) {
    using constants = aerobus::arithmetic_helpers<float>;
    float values[] = {
        0.078371367725F,
        0.707031250000F,
        0.708984375000F,
        constants::pi() / 8,
        3 * constants::pi() / 8,
        3 * constants::pi() / 4,
        3 * constants::pi() / 2,
        constants::pi(),
        9 * constants::pi() / 8,
        10.0F,
        -1.0F,
        776.0F,
    };

    for (float x : values) {
        float aero = aerobus::libm::cos(x);
        float expected = std::cos(x);
        EXPECT_TRUE((std::fabs(expected - aero) < std::numeric_limits<float>::epsilon())) << "aerobus::cos(" << x << ")"
            << " computed : " << std::hexfloat << aero << " but should " << expected << std::endl
            << "difference is : " << std::fabs(aero - expected) << std::endl;
    }

    float exact_values[14] = {
        std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN(),
        1E-26F, 1.0F,
        -1E-26F, 1.0F,
        -0.0F, 1.0F,
        0.0F, 1.0F,
        std::numeric_limits<float>::infinity(), std::numeric_limits<double>::quiet_NaN(),
        -std::numeric_limits<float>::infinity(), std::numeric_limits<double>::quiet_NaN()
    };

    for (unsigned long i = 0; i < sizeof(exact_values) / sizeof(float); i += 2) {
        float x = exact_values[i];
        float aero = aerobus::libm::cos(x);
        float expected = exact_values[i+1];
        if (x != x || expected != expected) {
            EXPECT_TRUE(aero != aero);
        } else {
            EXPECT_EQ(aero, expected) << "input was : " << x << std::endl;
        }
    }
}

int main() {
    return test_sin_exhaustive();
}