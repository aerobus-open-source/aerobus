namespace aerobus {
    struct i32 {
		using inner_type = int32_t;
		template<int32_t x>
		struct val {
			static constexpr int32_t v = x;

			template<typename valueType>
			static constexpr valueType get() { return static_cast<valueType>(x); }

			using is_zero_t = std::bool_constant<x == 0>;
			static std::string to_string() {
				return std::to_string(x);
			}

			template<typename valueRing>
			static constexpr valueRing eval(const valueRing& v) {
				return static_cast<valueRing>(x);
			}
		};

		using zero = val<0>;
		using one = val<1>;
		static constexpr bool is_field = false;
		static constexpr bool is_integral_domain = true;
		template<auto x>
		using inject_constant_t = val<static_cast<int32_t>(x)>;
		template<typename v>
		using inject_ring_t = v;

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
		struct remainder {
			using type = val<v1::v% v2::v>;
		};

		template<typename v1, typename v2>
		struct gt {
			using type = std::conditional_t<(v1::v > v2::v), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct lt {
			using type = std::conditional_t<(v1::v < v2::v), std::true_type, std::false_type>;
		};

		template<typename v1, typename v2>
		struct eq {
			using type = std::conditional_t<(v1::v == v2::v), std::true_type, std::false_type>;
		};

		template<typename v1>
		struct pos {
			using type = std::bool_constant<(v1::v > 0)>;
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
		using mod_t = typename remainder<v1, v2>::type;

		template<typename v1, typename v2>
		using gt_t = typename gt<v1, v2>::type;

		template<typename v1, typename v2>
		using lt_t = typename lt<v1, v2>::type;

		template<typename v1, typename v2>
		using eq_t = typename eq<v1, v2>::type;

		template<typename v1, typename v2>
		using gcd_t = typename internal::gcd<i32>::template type<v1, v2>;

		template<typename v>
		using pos_t = typename pos<v>::type;
	};
}