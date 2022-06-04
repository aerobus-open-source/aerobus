#include <cstdint>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>
#include <algorithm>

template <int64_t i, typename T, typename... Ts>
struct type_at
{
	static_assert(i < sizeof...(Ts) + 1, "index out of range");
	using type = typename type_at<i - 1, Ts...>::type;
};

template <typename T, typename... Ts> struct type_at<0, T, Ts...> {
	using type = T;
};

template <size_t i, typename... Ts>
using type_at_t = typename type_at<i, Ts...>::type;

template <std::size_t ... Is>
constexpr auto index_sequence_reverse(std::index_sequence<Is...> const&)
-> decltype(std::index_sequence<sizeof...(Is) - 1U - Is...>{});

template <std::size_t N>
using make_index_sequence_reverse
= decltype(index_sequence_reverse(std::make_index_sequence<N>{}));

template<typename Ring>
struct gcd {
	template<typename a, typename b, typename E = void>
	struct gcd_helper {};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<
			typename Ring::template eq_t<a, b>,
			std::true_type
		>::value
	>
	> {
		using type = a;
	};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<
			typename Ring::template gt_t<a, b>,
			std::true_type
		>::value
	>
	> {
		using type = typename gcd_helper<typename Ring::template sub_t<a, b>, b>::type;
	};

	template<typename a, typename b>
	struct gcd_helper<a, b, std::enable_if_t<
		std::is_same<
			typename Ring::template lt_t<a, b>,
			std::true_type
		>::value
	>
	> {
		using type = typename gcd_helper<a, typename Ring::template sub_t<b, a>>::type;
	};

	template<typename a, typename b>
	using type = typename gcd<Ring>::template gcd_helper<a, b>::type;
};

struct i32 {
	template<int32_t x>
	struct val {
		static constexpr int32_t v = x;
	};

	using zero = val<0>;
	using one = val<1>;

private:
	template<typename v1, typename v2>
	struct add {
		using type = val<v1::v + v2::v>;
	};

	template<typename v1, typename v2>
	struct sub {
		using type = val<v1::v - v2::v>;
	};

	template<typename v1, typename v2>
	struct mul {
		using type = val<v1::v* v2::v>;
	};

	template<typename v1, typename v2>
	struct div {
		using type = val<v1::v / v2::v>;
	};

	template<typename v1, typename v2>
	struct gt {
		using type = typename std::conditional<(v1::v > v2::v), std::true_type, std::false_type>::type;
	};

	template<typename v1, typename v2>
	struct lt {
		using type = typename std::conditional<(v1::v < v2::v), std::true_type, std::false_type>::type;
	};

	template<typename v1, typename v2>
	struct eq {
		using type = typename std::conditional<(v1::v == v2::v), std::true_type, std::false_type>::type;
	};

public:
	template<typename v1, typename v2>
	using add_t = typename add<v1, v2>::type;

	template<typename v1, typename v2>
	using sub_t = typename sub<v1, v2>::type;

	template<typename v1, typename v2>
	using mul_t = typename mul<v1, v2>::type;

	template<typename v1, typename v2>
	using div_t = typename div<v1, v2>::type;

	template<typename v1, typename v2>
	using gt_t = typename gt<v1, v2>::type;

	template<typename v1, typename v2>
	using lt_t = typename lt<v1, v2>::type;

	template<typename v1, typename v2>
	using eq_t = typename eq<v1, v2>::type;

	template<typename v1, typename v2>
	using gcd_t = typename gcd<i32>::template type<v1, v2>;
};

template<typename Ring, typename P, typename E = void>
struct poly_simplify;

template<typename Ring, typename P>
using poly_simplify_t = typename poly_simplify<Ring, P>::type;

template <typename Ring, typename P1, typename P2, typename I>
struct poly_add_low;

template <typename Ring, typename P1, typename P2, typename I>
struct poly_sub_low;

template<typename Ring, typename P1, typename P2>
using poly_add_t = poly_simplify_t<Ring, typename poly_add_low<
	Ring,
	P1,
	P2,
	make_index_sequence_reverse<
	std::max(P1::degree, P2::degree) + 1
	>>::type>;

template<typename Ring, typename P1, typename P2>
using poly_sub_t = poly_simplify_t<Ring, typename poly_sub_low<
	Ring,
	P1,
	P2,
	make_index_sequence_reverse<
	std::max(P1::degree, P2::degree) + 1
	>>::type>;

// coeffN x^N + ...
template<typename Ring>
struct polynomial {
	template<typename coeffN, typename... coeffs>
	struct val {
		static constexpr size_t degree = sizeof...(coeffs);
		using aN = coeffN;
		using strip = polynomial<Ring>::template val<coeffs...>;

		template<size_t index, typename E = void>
		struct coeff_at {};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index >= 0 && index <= sizeof...(coeffs))>> {
			using type = type_at_t<sizeof...(coeffs) - index, coeffN, coeffs...>;
		};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index < 0 || index > sizeof...(coeffs))>> {
			using type = typename Ring::zero;
		};

		template<size_t index>
		using coeff_at_t = typename coeff_at<index>::type;
	};

	// specialization for constants
	template<typename coeffN>
	struct val<coeffN> {
		static constexpr size_t degree = 0;
		using aN = coeffN;
		using strip = val<coeffN>;

		template<size_t index, typename E = void>
		struct coeff_at {};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index == 0)>> {
			using type = aN;
		};

		template<size_t index>
		struct coeff_at<index, std::enable_if_t<(index < 0 || index > 0)>> {
			using type = typename Ring::zero;
		};

		template<size_t index>
		using coeff_at_t = typename coeff_at<index>::type;
	};

private:
	template<typename v1, typename v2, typename E = void>
	struct eq_helper {};

	template<typename v1, typename v2>
	struct eq_helper<v1, v2, std::enable_if_t<v1::degree != v2::degree>> {
		using type = std::false_type;
	};

	template<typename v1, typename v2>
	struct eq_helper<v1, v2, std::enable_if_t<
								v1::degree == v2::degree && 
								v1::degree != 0 && 
								std::is_same<
									typename Ring::template eq_t<typename v1::aN, typename v2::aN>, 
									std::true_type
								>::value
							  >
	               > {
		using type = typename eq_helper<typename v1::strip, typename v2::strip>::type;
	};

	template<typename v1, typename v2>
	struct eq_helper<v1, v2, std::enable_if_t<
								v1::degree == v1::degree && 
								v1::degree == 0 
							 >
					> {
		using type = typename Ring::template eq_t<typename v1::aN, typename v2::aN>;
	};

public:
	using zero = typename polynomial<Ring>::template val<typename Ring::zero>;
	using one = typename polynomial<Ring>::template val<typename Ring::one>;

	template<typename v1, typename v2>
	using add_t = poly_add_t<Ring, v1, v2>;

	template<typename v1, typename v2>
	using sub_t = poly_sub_t<Ring, v1, v2>;

	template<typename v1, typename v2>
	using eq_t = typename eq_helper<v1, v2>::type;
};

template<typename Ring, typename P, typename E>
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

// addition at
template<typename Ring, typename P1, typename P2, size_t index>
struct poly_add_at {
	using type =
		typename Ring::template add_t<P1::template coeff_at_t<index>, P2::template coeff_at_t<index>>;
};

template<typename Ring, typename P1, typename P2, size_t index>
using poly_add_at_t = typename poly_add_at<Ring, P1, P2, index>::type;

template<typename Ring, typename P1, typename P2, std::size_t... I>
struct poly_add_low<Ring, P1, P2, std::index_sequence<I...>> {
	using type = typename polynomial<Ring>::template val<poly_add_at_t<Ring, P1, P2, I>...>;
};

// substraction at
template<typename Ring, typename P1, typename P2, size_t index>
struct poly_sub_at {
	using type =
		typename Ring::template sub_t<P1::template coeff_at_t<index>, P2::template coeff_at_t<index>>;
};

template<typename Ring, typename P1, typename P2, size_t index>
using poly_sub_at_t = typename poly_sub_at<Ring, P1, P2, index>::type;

template<typename Ring, typename P1, typename P2, std::size_t... I>
struct poly_sub_low<Ring, P1, P2, std::index_sequence<I...>> {
	using type = typename polynomial<Ring>::template val<poly_sub_at_t<Ring, P1, P2, I>...>;
};

template<typename Ring>
struct FractionField {
	template<typename val1, typename val2>
	struct val {
		using x = val1;
		using y = val2;
	};

	using zero = val<typename Ring::zero, typename Ring::one>;
	using one = val<typename Ring::one, typename Ring::one>;

	// for operators, we need gcd
};
