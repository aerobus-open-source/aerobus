
#include <string>
#include <type_traits>

namespace aerobus {
    namespace aad {
        template <typename derived_t, typename E = void> 
        struct base_traits;

        template<typename Derived>
        struct Expression {
            static std::string to_string() {
                return Derived::_to_string();
            }

            using value_at_t0 = typename base_traits<Derived>::value_at_t0;
            using derive_t = typename base_traits<Derived>::derive_t;
            using derive_x = typename base_traits<Derived>::derive_x;
        };
    
        // val is a Q64 type
        template<typename val>
        struct ConstantExpression: Expression<ConstantExpression<val>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "(" + val::to_string() + ")";
                }
        };

        // a*X^k
        template<typename a, size_t k>
        struct XExpression: Expression<XExpression<a, k>> {
            template<typename>
            friend struct Expression;

            private:
                static std::string _to_string() {
                    return "(" + a::to_string() + "X^" + std::to_string(k) + ")";
                }
        };

        // a*t^k
        template<typename a, size_t k>
        struct TExpression: Expression<TExpression<a, k>> {
            template<typename>
            friend struct Expression;

            private:
                static std::string _to_string() {
                    return "(" + a::to_string() + "t^" + std::to_string(k) + ")";
                }
        };

        template<typename E1, typename E2>
        struct AddExpression : Expression<AddExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "(" + E1::to_string() + " + " + E2::to_string() + ")";
                }
        };

        template<typename E1, typename E2>
        struct SubExpression : Expression<SubExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "(" + E1::to_string() + " - " + E2::to_string() + ")";
                }
        };
        
        template<typename E1, typename E2>
        struct MulExpression : Expression<MulExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "(" + E1::to_string() + " * " + E2::to_string() + ")";
                }
        };

        template<typename E1, typename E2>
        struct DivExpression : Expression<DivExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "(" + E1::to_string() + " / " + E2::to_string() + ")";
                }
        };

        // can't exist by itself but by composition
        template<typename E>
        struct ExpExpression : Expression<ExpExpression<E>> {
            template<typename>
            friend struct Expression;
            private:
                static std::string _to_string() {
                    return "exp(" + E::to_string() + ")";
                }
        };

        // template<typename E1, typename E2>
        // struct ComposeExpression: Expression<ComposeExpression<E1, E2>> {
        //     template<typename>
        //     friend struct Expression;
        //     private:
        //         static std::string _to_string() {
        //             return E1::to_string() + "(" + E2::to_string() + ")";
        //         }
        //     public:
        //         using _derive_t = MulExpression<
        //             ComposeExpression<typename E1::derive_t, E2>,
        //             typename E2::derive_t
        //         >;
        //         using _derive_x = MulExpression<
        //             ComposeExpression<typename E1::derive_x, E2>,
        //             typename E2::derive_t
        //         >;
        // };


        // base_traits specialization and simplification rules
        template<typename val>
        struct base_traits<ConstantExpression<val>> {
            using value_at_t0 = typename FPQ64::inject_ring_t<val>;
            using derive_x = ConstantExpression<val>;
            using derive_t = ConstantExpression<typename Q64::zero>;
        };

        template<typename a, size_t k>
        struct base_traits<XExpression<a, k>> {
            public:
                using value_at_t0 =  FPQ64::template val<polynomial<Q64>::monomial_t<a, k>, typename polynomial<Q64>::one>;
                using derive_t = ConstantExpression<typename Q64::zero>;
                using derive_x = std::conditional_t<(k == 0), ConstantExpression<typename Q64::zero>,
                                    std::conditional_t<k == 1, ConstantExpression<a>, 
                                        XExpression<Q64::mul_t<a, Q64::inject_constant_t<k>>, k-1>>>;
        };
        
        template<typename a, size_t k>
        struct base_traits<TExpression<a, k>> {
            using value_at_t0 = std::conditional_t<
                k == 0,
                typename FPQ64::inject_ring_t<a>,
                typename FPQ64::zero
            >;
            using derive_t = std::conditional_t<(k == 0), 
                                ConstantExpression<typename Q64::zero>,
                                std::conditional_t<k == 1, 
                                        ConstantExpression<a>, 
                                        TExpression<Q64::mul_t<a, Q64::inject_constant_t<k>>, k-1>>>;
            using derive_x = ConstantExpression<typename Q64::zero>;
        };

        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>> {
            using value_at_t0 = FPQ64::template add_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = AddExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = AddExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>> {
            using value_at_t0 = FPQ64::template sub_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = SubExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = SubExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>> {
            using value_at_t0 = FPQ64::template mul_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = AddExpression<
                            MulExpression<typename base_traits<E1>::derive_t, E2>,
                            MulExpression<E1, typename base_traits<E2>::derive_t>>;
            
            using derive_x = AddExpression<
                            MulExpression<typename base_traits<E1>::derive_x, E2>,
                            MulExpression<E1, typename base_traits<E2>::derive_x>>;
        };

        template<typename E1, typename E2>
        struct base_traits<DivExpression<E1, E2>> {
            using value_at_t0 = typename FPQ64::template div_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = DivExpression<
                                SubExpression<
                                    MulExpression<typename base_traits<E1>::derive_t, E2>,
                                    MulExpression<E1, typename base_traits<E2>::derive_t>
                                >,
                                MulExpression<E2, E2>
                            >;
            
            using derive_x = DivExpression<
                                SubExpression<
                                    MulExpression<typename base_traits<E1>::derive_x, E2>,
                                    MulExpression<E1, typename base_traits<E2>::derive_x>
                                >,
                                MulExpression<E2, E2>
                            >;
        };

        template<typename E>
        struct base_traits<ExpExpression<E>> {
            using value_at_t0 = typename FPQ64::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, ExpExpression<E>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, ExpExpression<E>>;
        };

        template<typename E, size_t k>
        struct taylor_coeff_helper {
            using type = FPQ64::template div_t<typename taylor_coeff_helper<typename base_traits<E>::derive_t, k-1>::type,
                                                FPQ64::inject_constant_t<k>> ;
        };

        template<typename E>
        struct taylor_coeff_helper<E, 0> {
            using type = typename base_traits<E>::value_at_t0;
        };

        template<typename E, size_t k>
        using taylor_coeff_t = typename taylor_coeff_helper<E, k>::type;

        template<typename E, size_t k>
        struct taylor_expansion_helper {
            using type = FPQ64::template add_t<taylor_coeff_t<E, k>, typename taylor_expansion_helper<E, k-1>::type>;
        };

        template<typename E>
        struct taylor_expansion_helper<E, 0> {
            using type = taylor_coeff_t<E, 0>;
        };

        template<typename E, size_t k>
        using taylor_expansion_t = typename taylor_expansion_helper<E, k>::type;
    }
}