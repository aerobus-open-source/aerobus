#include <gtest/gtest.h>
#include <cstdint>

#include "aerobus.h"

TEST(memory, aligned_malloc_and_free) {
    // 3 doubles = 24 bytes, not a multiple of the 64-byte alignment
    double* p = aerobus::aligned_malloc<double>(3, 64);
    ASSERT_NE(p, nullptr);
    EXPECT_EQ(reinterpret_cast<std::uintptr_t>(p) % 64, 0u);
    p[0] = 1.0;
    p[2] = 3.0;
    EXPECT_EQ(p[0] + p[2], 4.0);
    aerobus::aligned_free(p);
}

TEST(memory, aligned_malloc_size_overflow) {
    // count * sizeof(double) = 2^64 + 8 wraps to 8: must fail, not return a 64-byte block (#49)
    const size_t n = (SIZE_MAX >> 3) + 2;
    double* p = aerobus::aligned_malloc<double>(n, 64);
    EXPECT_EQ(p, nullptr);
    aerobus::aligned_free(p);
    // rounding count * sizeof(T) up to the alignment must not wrap either
    char* q = aerobus::aligned_malloc<char>(SIZE_MAX - 10, 64);
    EXPECT_EQ(q, nullptr);
    aerobus::aligned_free(q);
    EXPECT_EQ(aerobus::aligned_malloc<double>(1, 0), nullptr);
}
