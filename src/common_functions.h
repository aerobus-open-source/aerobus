namespace aerobus {
    #define MAKE_Q64(x) Q64::inject_constant_t<x>

	template<typename T, size_t x, typename E = void>
	struct factorial {};

	template<typename T, size_t x>
	struct factorial<T, x, std::enable_if_t<(x > 0)>> {
	private:
		template<typename, size_t, typename>
		friend struct factorial;
	public:
		using type = typename T::template mul_t<typename T::template val<x>, typename factorial<T, x - 1>::type>;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T>
	struct factorial<T, 0> {
	public:
		using type = typename T::one;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	namespace internal {
		template<typename T, size_t k, size_t n, typename E = void>
		struct combination_helper {};

		template<typename T, size_t k, size_t n>
		struct combination_helper<T, k, n, std::enable_if_t<(n >= 0 && k <= (n / 2) && k > 0)>> {
			using type = typename FractionField<T>::template mul_t<
				typename combination_helper<T, k - 1, n - 1>::type,
				typename FractionField<T>::template val<typename T::template val<n>, typename T::template val<k>>>;
		};

		template<typename T, size_t k, size_t n>
		struct combination_helper<T, k, n, std::enable_if_t<(n >= 0 && k > (n / 2) && k > 0)>> {
			using type = typename combination_helper<T, n - k, n>::type;
		};

		template<typename T, size_t n>
		struct combination_helper<T, 0, n> {
			using type = typename FractionField<T>::one;
		};
	}

	template<typename T, size_t k, size_t n>
	struct combination {
		using type = typename internal::combination_helper<T, k, n>::type::x;
		static constexpr typename T::inner_type value = internal::combination_helper<T, k, n>::type::template get<typename T::inner_type>();
	};


	template<typename T, size_t m>
	struct bernouilli;

	namespace internal {
		template<typename T, typename accum, size_t k, size_t m>
		struct bernouilli_helper {
			using type = typename bernouilli_helper<
				T,
				typename FractionField<T>::template add_t<
				accum,
				typename FractionField<T>::template mul_t<
				typename FractionField<T>::template val<
				typename combination<T, k, m + 1>::type,
				typename T::one>,
				typename bernouilli<T, k>::type
				>
				>,
				k + 1,
				m
			>::type;
		};

		template<typename T, typename accum, size_t m>
		struct bernouilli_helper<T, accum, m, m>
		{
			using type = accum;
		};
	}

	template<typename T, size_t m>
	struct bernouilli {
		using type = typename FractionField<T>::template mul_t<
			typename internal::bernouilli_helper<T, typename FractionField<T>::zero, 0, m>::type,
			typename FractionField<T>::template val<
			typename T::template val<static_cast<typename T::inner_type>(-1)>,
			typename T::template val<static_cast<typename T::inner_type>(m + 1)>
			>
		>;

		template<typename floatType>
		static constexpr floatType value = type::template get<floatType>();
	};

	template<typename T>
	struct bernouilli<T, 0> {
		using type = typename FractionField<T>::one;

		template<typename floatType>
		static constexpr floatType value = type::template get<floatType>();
	};

	template<typename T, int k, typename E = void>
	struct alternate {};

	template<typename T, int k>
	struct alternate<T, k, std::enable_if_t<k % 2 == 0>> {
		using type = typename T::one;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T, int k>
	struct alternate<T, k, std::enable_if_t<k % 2 != 0>> {
		using type = typename T::template sub_t<typename T::zero, typename T::one>;
		static constexpr typename T::inner_type value = type::template get<typename T::inner_type>();
	};

	template<typename T, auto p, auto n>
	struct pow {
		using type = typename T::template mul_t<typename T::template val<p>, typename pow<T, p, n - 1>::type>;
	};

	template<typename T, auto p>
	struct pow<T, p, 0> { using type = typename T::one; };

	namespace internal {
		template<typename, template<typename, size_t> typename, class>
		struct make_taylor_impl;

		template<typename T, template<typename, size_t> typename coeff_at, size_t... Is>
		struct make_taylor_impl<T, coeff_at, std::integer_sequence<size_t, Is...>> {
			using type = typename polynomial<FractionField<T>>::template val<typename coeff_at<T, Is>::type...>;
		};
	}

	// generic taylor serie, depending on coefficients
	template<typename T, template<typename, size_t index> typename coeff_at, size_t deg>
	using taylor = typename internal::make_taylor_impl<T, coeff_at, internal::make_index_sequence_reverse<deg + 1>>::type;

	namespace internal {
		template<typename T, size_t i>
		struct exp_coeff {
			using type = typename FractionField<T>::template val<typename T::one, typename factorial<T, i>::type>;
		};

		template<typename T, size_t i, typename E = void>
		struct sin_coeff_helper {};

		template<typename T, size_t i>
		struct sin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct sin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
			using type = typename FractionField<T>::template val<typename alternate<T, i / 2>::type, typename factorial<T, i>::type>;
		};

		template<typename T, size_t i>
		struct sin_coeff {
			using type = typename sin_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i, typename E = void>
		struct sh_coeff_helper {};

		template<typename T, size_t i>
		struct sh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct sh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
			using type = typename FractionField<T>::template val<typename T::one, typename factorial<T, i>::type>;
		};

		template<typename T, size_t i>
		struct sh_coeff {
			using type = typename sh_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i, typename E = void>
		struct cos_coeff_helper {};

		template<typename T, size_t i>
		struct cos_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct cos_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
			using type = typename FractionField<T>::template val<typename alternate<T, i / 2>::type, typename factorial<T, i>::type>;
		};

		template<typename T, size_t i>
		struct cos_coeff {
			using type = typename cos_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i, typename E = void>
		struct cosh_coeff_helper {};

		template<typename T, size_t i>
		struct cosh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct cosh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
			using type = typename FractionField<T>::template val<typename T::one, typename factorial<T, i>::type>;
		};

		template<typename T, size_t i>
		struct cosh_coeff {
			using type = typename cosh_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i>
		struct geom_coeff { using type = typename FractionField<T>::one; };


		template<typename T, size_t i, typename E = void>
		struct atan_coeff_helper;

		template<typename T, size_t i>
		struct atan_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type> {
			using type = typename FractionField<T>::template val<alternate<T, i / 2>::type, typename T::template val<i>>;
		};

		template<typename T, size_t i>
		struct atan_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct atan_coeff { using type = typename atan_coeff_helper<T, i>::type; };

		template<typename T, size_t i, typename E = void>
		struct asin_coeff_helper;

		template<typename T, size_t i>
		struct asin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type>
		{
			using type = typename FractionField<T>::template val<
				typename factorial<T, i - 1>::type,
				typename T::template mul_t<
				typename T::template val<i>,
				typename T::template mul_t<
				typename pow<T, 4, i / 2>::type,
				pow<T, factorial<T, i / 2>::value, 2>
				>>>;
		};

		template<typename T, size_t i>
		struct asin_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type>
		{
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct asin_coeff {
			using type = typename asin_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i>
		struct lnp1_coeff {
			using type = typename FractionField<T>::template val<
				typename alternate<T, i + 1>::type,
				typename T::template val<i>>;
		};

		template<typename T>
		struct lnp1_coeff<T, 0> { using type = typename FractionField<T>::zero; };

		template<typename T, size_t i, typename E = void>
		struct asinh_coeff_helper;

		template<typename T, size_t i>
		struct asinh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type>
		{
			using type = typename FractionField<T>::template val<
				typename T::template mul_t<
				typename alternate<T, i / 2>::type,
				typename factorial<T, i - 1>::type
				>,
				typename T::template mul_t<
				typename T::template mul_t<
				typename T::template val<i>,
				typename pow<T, (factorial<T, i / 2>::value), 2>::type
				>,
				typename pow<T, 4, i / 2>::type
				>
			>;
		};

		template<typename T, size_t i>
		struct asinh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type>
		{
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct asinh_coeff {
			using type = typename asinh_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i, typename E = void>
		struct atanh_coeff_helper;

		template<typename T, size_t i>
		struct atanh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 1>::type>
		{
			// 1/i
			using type = typename FractionField<T>:: template val<
				typename T::one,
				typename T::template val<static_cast<typename T::inner_type>(i)>>;
		};

		template<typename T, size_t i>
		struct atanh_coeff_helper<T, i, typename std::enable_if<(i & 1) == 0>::type>
		{
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct atanh_coeff {
			using type = typename asinh_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i, typename E = void>
		struct tan_coeff_helper;

		template<typename T, size_t i>
		struct tan_coeff_helper<T, i, typename std::enable_if_t<(i % 2) == 0>> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct tan_coeff_helper<T, i, typename std::enable_if_t<(i % 2) != 0>> {
		private:
			// 4^((i+1)/2)
			using _4p = typename FractionField<T>::template inject_t<typename pow<T, 4, (i + 1) / 2>::type>;
			// 4^((i+1)/2) - 1
			using _4pm1 = typename FractionField<T>::template sub_t<_4p, typename FractionField<T>::one>;
			// (-1)^((i-1)/2)
			using altp = typename FractionField<T>::template inject_t<typename alternate<T, (i - 1) / 2>::type>;
			using dividend = typename FractionField<T>::template mul_t<
				altp,
				typename FractionField<T>::template mul_t<
				_4p,
				typename FractionField<T>::template mul_t<
				_4pm1,
				typename bernouilli<T, (i + 1)>::type
				>
				>
			>;
		public:
			using type = typename FractionField<T>::template div_t<dividend,
				typename FractionField<T>::template inject_t<typename factorial<T, i + 1>::type>>;
		};

		template<typename T, size_t i>
		struct tan_coeff {
			using type = typename tan_coeff_helper<T, i>::type;
		};

		template<typename T, size_t i, typename E = void>
		struct tanh_coeff_helper;

		template<typename T, size_t i>
		struct tanh_coeff_helper<T, i, typename std::enable_if_t<(i % 2) == 0>> {
			using type = typename FractionField<T>::zero;
		};

		template<typename T, size_t i>
		struct tanh_coeff_helper<T, i, typename std::enable_if_t<(i % 2) != 0>> {
		private:
			using _4p = typename FractionField<T>::template inject_t<typename pow<T, 4, (i + 1) / 2>::type>;
			using _4pm1 = typename FractionField<T>::template sub_t<_4p, typename FractionField<T>::one>;
			using dividend =
				typename FractionField<T>::template mul_t<
				_4p,
				typename FractionField<T>::template mul_t<
				_4pm1,
				typename bernouilli<T, (i + 1)>::type
				>
				>::type;
		public:
			using type = typename FractionField<T>::template div_t<dividend,
				typename FractionField<T>::template inject_t<typename factorial<T, i + 1>::type>>;
		};

		template<typename T, size_t i>
		struct tanh_coeff {
			using type = typename tanh_coeff_helper<T, i>::type;
		};
	}

	// e^x
	template<typename T, size_t deg>
	using exp = taylor<T, internal::exp_coeff, deg>;

	// e^x - 1
	template<typename T, size_t deg>
	using expm1 = typename polynomial<FractionField<T>>::template sub_t<
		exp<T, deg>,
		typename polynomial<FractionField<T>>::one>;

	/// ln(1+x)
	template<typename T, size_t deg>
	using lnp1 = taylor<T, internal::lnp1_coeff, deg>;

	/// atan(x);
	template<typename T, size_t deg>
	using atan = taylor<T, internal::atan_coeff, deg>;

	/// sin(x)
	template<typename T, size_t deg>
	using sin = taylor<T, internal::sin_coeff, deg>;

	/// sh(x)
	template<typename T, size_t deg>
	using sinh = taylor<T, internal::sh_coeff, deg>;

	/// ch(x)
	template<typename T, size_t deg>
	using cosh = taylor<T, internal::cosh_coeff, deg>;

	/// cos(x)
	template<typename T, size_t deg>
	using cos = taylor<T, internal::cos_coeff, deg>;

	/// 1 / (1-x)
	template<typename T, size_t deg>
	using geometric_sum = taylor<T, internal::geom_coeff, deg>;

	/// asin(x)
	template<typename T, size_t deg>
	using asin = taylor<T, internal::asin_coeff, deg>;

	/// asinh(x)
	template<typename T, size_t deg>
	using asinh = taylor<T, internal::asinh_coeff, deg>;

	/// atanh(x)
	template<typename T, size_t deg>
	using atanh = taylor<T, internal::atanh_coeff, deg>;

	/// tan(x)
	template<typename T, size_t deg>
	using tan = taylor<T, internal::tan_coeff, deg>;

	/// tanh(x)
	template<typename T, size_t deg>
	using tanh = taylor<T, internal::tanh_coeff, deg>;

}