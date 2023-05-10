
#ifdef _MSC_VER
#define ALIGNED(x) __declspec(align(x))
#define INLINED __forceinline 
#else 
#define ALIGNED(x) __attribute__((aligned(x)))
#define INLINED __attribute__((always_inline))  
#endif
namespace aerobus {

	namespace internal
	{
		template<template<typename...> typename TT, typename T>
		struct is_instantiation_of : std::false_type { };

		template<template<typename...> typename TT, typename... Ts>
		struct is_instantiation_of<TT, TT<Ts...>> : std::true_type { };

		template<template<typename...> typename TT, typename T>
		inline constexpr bool is_instantiation_of_v = is_instantiation_of<TT, T>::value;

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


		template<int32_t n, int32_t i, typename E = void>
		struct _is_prime {};

		template<int32_t n, int32_t i>
		struct _is_prime<n, i, std::enable_if_t<(n == 2 || n == 3)>> : std::true_type {};

		template<int32_t n, int32_t i>
		struct _is_prime<n, i, std::enable_if_t<(n != 2 && n != 3) && (n <= 1 || n % 2 == 0 || n % 3 == 0)>> : std::false_type {};

		template<int32_t n, int32_t i>
		struct _is_prime<n, i, std::enable_if_t<
			(n > 1 && n % 2 != 0 && n % 3 != 0) &&
			(i >= 5 && i * i <= n) &&
			(n% i == 0 || n % (i + 2) == 0)>> : std::false_type {};


		template<int32_t n, int32_t i>
		struct _is_prime<n, i, std::enable_if_t<
			(n > 1 && n % 2 != 0 && n % 3 != 0) &&
			(i >= 5 && i * i <= n) &&
			(n% i != 0 && n % (i + 2) != 0)>> {
			static constexpr bool value = _is_prime<n, i + 6>::value;
		};

		template<int32_t n, int32_t i>
		struct _is_prime<n, i, std::enable_if_t<
			(n > 1 && n % 2 != 0 && n % 3 != 0) &&
			(i >= 5 && i * i > n)>> : std::true_type {};
	}

	template<int32_t n>
	struct is_prime {
		static constexpr bool value = internal::_is_prime<n, 5>::value;
	};

	namespace internal {
		template <std::size_t... Is>
		constexpr auto index_sequence_reverse(std::index_sequence<Is...> const&)
			-> decltype(std::index_sequence<sizeof...(Is) - 1U - Is...>{});

		template <std::size_t N>
		using make_index_sequence_reverse
			= decltype(index_sequence_reverse(std::make_index_sequence<N>{}));

		template<typename Ring, typename E = void>
		struct gcd;

		template<typename Ring>
		struct gcd<Ring, std::enable_if_t<Ring::is_integral_domain>> {
			template<typename A, typename B, typename E = void>
			struct gcd_helper {};

			// B = 0, A > 0
			template<typename A, typename B>
			struct gcd_helper<A, B, std::enable_if_t<
				((B::is_zero_t::value) &&
					(Ring::template gt_t<A, typename Ring::zero>::value))>>
			{
				using type = A;
			};

			// B = 0, A < 0
			template<typename A, typename B>
			struct gcd_helper<A, B, std::enable_if_t<
				((B::is_zero_t::value) &&
					!(Ring::template gt_t<A, typename Ring::zero>::value))>>
			{
				using type = typename Ring::template sub_t<typename Ring::zero, A>;
			};

			// B != 0
			template<typename A, typename B>
			struct gcd_helper<A, B, std::enable_if_t<
				(!B::is_zero_t::value)
				>> {
			private:
				// A / B
				using k = typename Ring::template div_t<A, B>;
				// A - (A/B)*B = A % B
				using m = typename Ring::template sub_t<A, typename Ring::template mul_t<k, B>>;
			public:
				using type = typename gcd_helper<B, m>::type;
			};

			template<typename A, typename B>
			using type = typename gcd_helper<A, B>::type;
		};
	}
}
