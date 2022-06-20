#include <cstdint> // NOLINT(clang-diagnostic-pragma-pack)
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <functional>
#include <string>

#include "utilities.h"
#include "i32.h"
#include "i64.h"
#include "zpz.h"
#include "polynomial.h"
#include "fraction_field.h"
namespace aerobus {
	using Q32 = FractionField<i32>;
	using FPQ32 = FractionField<polynomial<Q32>>;
	using Q64 = FractionField<i64>;
	using FPQ64 = FractionField<polynomial<Q64>>;
}
#include "common_functions.h"
#include "aad.h"
#include "generating_functions.h"