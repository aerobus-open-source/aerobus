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
