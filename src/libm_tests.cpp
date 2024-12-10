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

int main() {
    return test_sin_exhaustive();
}