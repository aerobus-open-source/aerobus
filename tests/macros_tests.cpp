#include <gtest/gtest.h>

#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"

// aerobus.h must not leak macros with generic names into user code
#if defined(INLINED) || defined(DEVICE) || defined(ALIGNED) || defined(ZPZV) || defined(POLYV)
#error "aerobus.h leaks an unprefixed macro"
#endif

// nor use an identifier reserved to the implementation as include guard
#ifdef __INC_AEROBUS__
#error "aerobus.h uses a reserved include guard"
#endif

#if !defined(AEROBUS_INLINE) || !defined(AEROBUS_DEVICE) || !defined(AEROBUS_ALIGNED)
#error "aerobus.h must define the AEROBUS_ prefixed macros"
#endif

TEST(api, only_prefixed_macros) {
    SUCCEED();
}
