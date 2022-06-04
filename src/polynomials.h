#pragma once
#include <utility>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include "type_utils.h"

// coeffN x^N + ...
template<typename Ring>
struct polynomial {
	template<typename coeffN, typename... coeffs>
	struct val {
		static constexpr size_t degree = sizeof...(coeffs);
		using aN = coeffN;
		using strip = std::conditional_t<
			(sizeof...(coeffs) > 0), 
			typename polynomial<Ring>::template val<coeffs...>, 
			typename polynomial<Ring>::template val<coeffN, coeffs...>>;

		template<size_t index>
		using coeff_at = std::conditional_t<
			(index < 0 || index > sizeof...(coeffs)), 
			typename Ring::zero, 
			type_at_t<sizeof...(coeffs) - index, coeffN, coeffs...>>;
	};

	using zero = typename polynomial<Ring>::template val<typename Ring::zero>;
	using one = typename polynomial<Ring>::template val<typename Ring::one>;
};

template<typename Ring, typename P1, typename P2>
struct poly_add {
	using type = polynomial<Ring>::template val<typename Ring::zero>;
};

template<typename Ring, typename P, typename E = void>
struct poly_simplify {};

// when high power is zero : strip
template<typename Ring, typename P>
struct poly_simplify<Ring, P, typename std::enable_if<
	std::is_same<
	typename Ring::zero,
	typename P::aN
	>::value && (P::degree > 0)
>::type>
{
	using type = typename poly_simplify<Ring, typename P::strip>::type;
};

// otherwise : do nothing
template<typename Ring, typename P>
struct poly_simplify<Ring, P, typename std::enable_if<
	!std::is_same<
	typename Ring::zero,
	typename P::aN
	>::value && (P::degree > 0)
>::type>
{
	using type = P;
};

// do not simplify constants
template<typename Ring, typename P>
struct poly_simplify<Ring, P, std::enable_if_t<P::degree == 0>> {
	using type = P;
};


template<typename Ring, typename P>
using poly_simplify_t = poly_simplify<Ring, P>::type;
