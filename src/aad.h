#pragma once

namespace aerobus {
    using P64X = polynomial<Q64, 'x'>;
    using FPQ64X = RationalFraction<Q64, 'x'>;
    using P64T = polynomial<P64X, 't'>;

    namespace aad {
        template <typename derived_t, typename E = void> 
        struct base_traits;

        template<typename E>
        using simplify_t = typename base_traits<E>::type;

        template<typename E1, typename E2, typename E = void>
        struct eq_helper : std::false_type {};

        template<typename E1, typename E2>
        struct eq_helper<E1, E2, 
            std::enable_if_t<
                std::is_same<simplify_t<E1>, simplify_t<E2>>::value>>: std::true_type {};

        template<typename E1, typename E2>
        using eq_t = typename eq_helper<E1, E2>::type;

        template<typename Derived>
        struct Expression {
            using type = typename base_traits<Derived>::type;
            using value_at_t0 = typename base_traits<Derived>::value_at_t0;
            using derive_t = typename base_traits<Derived>::derive_t;
            using derive_x = typename base_traits<Derived>::derive_x;

            static std::string to_string() {
                return type::_to_string();
            }
        };
    
        // val is a Q64 type
        template<typename val>
        struct ConstantExpression: Expression<ConstantExpression<val>> {
            template<typename>
            friend struct Expression;
            using v = val;
            private:
                static std::string _to_string() {
                    return "(" + val::to_string() + ")";
                }
        };

        template<typename E>
        struct MinExpression: Expression<MinExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;

            private:
                static std::string _to_string() {
                    return "-(" + E::to_string() + ")";
                }
        };

        // a*X^k
        template<typename a, size_t k>
        struct XExpression: Expression<XExpression<a, k>> {
            template<typename>
            friend struct Expression;

            private:
                static std::string _to_string() {
                    std::string result = "";
                    if (!a::is_zero_t::value && !std::is_same<a, typename Q64::one>::value) {
                        result += a::to_string();
                    }
                    result += "X";
                    if (k > 1) {
                        result += "^" + std::to_string(k);
                    }
                    return result;
                }
        };

        // a*t^k
        template<typename a, size_t k>
        struct TExpression: Expression<TExpression<a, k>> {
            template<typename>
            friend struct Expression;

            private:
                static std::string _to_string() {
                    std::string result = "";
                    if (!a::is_zero_t::value && !std::is_same<a, typename Q64::one>::value) {
                        result += a::to_string();
                    }
                    result += "t";
                    if (k > 1) {
                        result += "^" + std::to_string(k);
                    }
                    return result;
                }
        };

        template<typename E1, typename E2>
        struct AddExpression : Expression<AddExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            using left = E1;
            using right = E2;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") + (" + E2::to_string() + "))";
                }
        };

        template<typename E1, typename E2>
        struct SubExpression : Expression<SubExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            using left = E1;
            using right = E2;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") - (" + E2::to_string() + "))";
                }
        };
        
        template<typename E1, typename E2>
        struct MulExpression : Expression<MulExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            using left = E1;
            using right = E2;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") * (" + E2::to_string() + "))";
                }
        };

        template<typename E1, typename E2>
        struct DivExpression : Expression<DivExpression<E1, E2>> {
            template<typename>
            friend struct Expression;
            using left = E1;
            using right = E2;
            private:
                static std::string _to_string() {
                    return "((" + E1::to_string() + ") / (" + E2::to_string() + "))";
                }
        };

        // can't exist by itself but by composition
        template<typename E>
        struct ExpExpression : Expression<ExpExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;
            private:
                static std::string _to_string() {
                    return "exp(" + E::to_string() + ")";
                }
        }; 
        
        template<typename E>
        struct SinExpression : Expression<SinExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;
            private:
                static std::string _to_string() {
                    return "sin(" + E::to_string() + ")";
                }
        };
        
        template<typename E>
        struct CosExpression : Expression<CosExpression<E>> {
            template<typename>
            friend struct Expression;
            using v = E;
            private:
                static std::string _to_string() {
                    return "cos(" + E::to_string() + ")";
                }
        };

        template<typename E1, typename E2>
        struct base_traits<DivExpression<E1, E2>> {
            using type = DivExpression<E1, E2>;
            using value_at_t0 = typename FPQ64X::template div_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
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

        template<typename val>
        struct base_traits<ConstantExpression<val>> {
            using type = ConstantExpression<val>;
            using value_at_t0 = typename FPQ64X::inject_ring_t<val>;
            using derive_x = ConstantExpression<val>;
            using derive_t = ConstantExpression<typename Q64::zero>;
        };

        template<typename a, size_t k>
        struct base_traits<XExpression<a, k>> {
            public:
                using type = XExpression<a, k>;
                using value_at_t0 =  FPQ64X::template val<polynomial<Q64>::monomial_t<a, k>, typename polynomial<Q64>::one>;
                using derive_t = ConstantExpression<typename Q64::zero>;
                using derive_x = std::conditional_t<(k == 0), ConstantExpression<typename Q64::zero>,
                                    std::conditional_t<k == 1, ConstantExpression<a>, 
                                        XExpression<Q64::mul_t<a, Q64::inject_constant_t<k>>, k-1>>>;
        };
        
        template<typename a, size_t k>
        struct base_traits<TExpression<a, k>> {
            using type = TExpression<a, k>;
            using value_at_t0 = std::conditional_t<
                k == 0,
                typename FPQ64X::inject_ring_t<a>,
                typename FPQ64X::zero
            >;
            using derive_t = std::conditional_t<(k == 0), 
                                ConstantExpression<typename Q64::zero>,
                                std::conditional_t<k == 1, 
                                        ConstantExpression<a>, 
                                        TExpression<Q64::mul_t<a, Q64::inject_constant_t<k>>, k-1>>>;
            using derive_x = ConstantExpression<typename Q64::zero>;
        };

        template<typename E>
        struct base_traits<ExpExpression<E>, 
            std::enable_if_t<!eq_t<ConstantExpression<typename Q64::zero>, simplify_t<E>>::value>
        > {
            using type = ExpExpression<E>;
            using value_at_t0 = typename FPQ64X::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, ExpExpression<simplify_t<E>>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, ExpExpression<simplify_t<E>>>;
        };
        
        template<typename E>
        struct base_traits<SinExpression<E>> {
            using type = SinExpression<E>;
            using value_at_t0 = typename FPQ64X::zero; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, CosExpression<E>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, CosExpression<E>>;
        };
        
        template<typename E>
        struct base_traits<CosExpression<E>> {
            using type = CosExpression<E>;
            using value_at_t0 = typename FPQ64X::one; // TODO: check that E::value_at_0 == 0
            using derive_t = MulExpression<typename base_traits<E>::derive_t, SubExpression<ConstantExpression<Q64::zero>, SinExpression<E>>>;
            using derive_x = MulExpression<typename base_traits<E>::derive_x, SubExpression<ConstantExpression<Q64::zero>, SinExpression<E>>>;
        };

        // export facilities
        using T = aad::TExpression<typename Q64::one, 1>;
        using X = aad::XExpression<typename  Q64::one, 1>;
        using XT = aad::MulExpression<X, T>;
        using ONE = aad::ConstantExpression<typename Q64::one>;
        using ZERO = aad::ConstantExpression<typename Q64::zero>;
        using TWO = aad::ConstantExpression<Q64::inject_constant_t<2>>;

        // equalities
        template<typename E1, typename E2>
        struct eq_helper<MulExpression<E1, E2>, MulExpression<E2, E1>> : std::true_type {};

        template<typename E1, typename E2>
        struct eq_helper<AddExpression<E1, E2>, AddExpression<E2, E1>> : std::true_type {};

        template<typename E1, typename E2, typename E3>
        struct eq_helper<AddExpression<AddExpression<E1, E2>, E3>, AddExpression<E1, AddExpression<E2, E3>>> : std::true_type {};

        template<typename E1, typename E2, typename E3>
        struct eq_helper<MulExpression<MulExpression<E1, E2>, E3>, MulExpression<E1, MulExpression<E2, E3>>> : std::true_type {};

        template<typename E1, typename E2, typename E3>
        struct eq_helper<
            MulExpression<E1, AddExpression<E2, E3>>, 
            AddExpression<
                MulExpression<E1, E2>, 
                MulExpression<E1, E3>
            >
        > : std::true_type {};

        
        template<typename E1, typename E2>
        struct eq_helper<
            SubExpression<ZERO, SubExpression<E1, E2>>, 
            SubExpression<E2, E1>
            > : std::true_type {};

        template<typename E1, typename E2, typename E3>
        struct eq_helper<
            MulExpression<AddExpression<E2, E3>, E1>, 
            AddExpression<
                MulExpression<E2, E1>, 
                MulExpression<E3, E1>
            >
        > : std::true_type {};

        template<typename E1, typename E2, typename E3>
        struct eq_helper<
            MulExpression<E1, SubExpression<E2, E3>>, 
            SubExpression<
                MulExpression<E1, E2>, 
                MulExpression<E1, E3>
            >
        > : std::true_type {};

        template<typename E1, typename E2, typename E3>
        struct eq_helper<
            MulExpression<SubExpression<E2, E3>, E1>, 
            SubExpression<
                MulExpression<E2, E1>, 
                MulExpression<E3, E1>
            >
        > : std::true_type {};

        template<typename E>
        struct eq_helper<MinExpression<MinExpression<E>>, E>: std::true_type {};

        template<typename E>
        struct eq_helper<E, MinExpression<MinExpression<E>>>: std::true_type {};

        template<template<typename E> typename expr, typename E1, typename E2>
        struct eq_helper<
                expr<E1>, expr<E2>,
                std::enable_if_t<!std::is_same<expr<E1>, expr<E2>>::value>
                >: eq_helper<E1, E2> {};

        struct simplification_rules {
            template <typename E1, typename E2>
            using zero_right =
                std::conjunction<
                    std::negation<internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>>,
                    eq_t<ZERO, simplify_t<E2>>>;

            template <typename E1, typename E2>
            using one_right =
                std::conjunction<
                    std::negation<internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>>,
                    eq_t<ONE, simplify_t<E2>>>;
            
            template <typename E1, typename E2>
            using zero_left = zero_right<E2, E1>;
            
            template <typename E1, typename E2>
            using one_left = one_right<E2, E1>;

            template<typename E1, typename E2>
            using right_min_left = std::conjunction<
                internal::is_instantiation_of<MinExpression, E2>,
                eq_t<E2, MinExpression<E1>>,
                std::negation<internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>>
            >;
            
            template<typename E1, typename E2>
            using left_min_right = std::conjunction<
                internal::is_instantiation_of<MinExpression, E1>,
                eq_t<E1, MinExpression<E2>>,
                std::negation<internal::is_instantiation_of<ConstantExpression, simplify_t<E2>>>
            >;

            template<typename E1, typename E2>
            using two_constants = 
                std::conjunction<
                    internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>,
                    internal::is_instantiation_of<ConstantExpression, simplify_t<E2>>
                >;

            template<typename E1, typename E2>
            using two_exp = 
                std::conjunction<
                    internal::is_instantiation_of<ExpExpression, simplify_t<E1>>,
                    internal::is_instantiation_of<ExpExpression, simplify_t<E2>>
                >;

            template<typename E1, typename E2>
            using sub_right = std::conjunction<
                internal::is_instantiation_of<SubExpression, simplify_t<E2>>,
                std::negation<one_left<E1, E2>>,
                std::negation<zero_left<E1, E2>>
            >; 
            
            template<typename E1, typename E2>
            using sub_left = std::conjunction<
                internal::is_instantiation_of<SubExpression, simplify_t<E1>>,
                std::negation<internal::is_instantiation_of<AddExpression, simplify_t<E2>>>,
                std::negation<internal::is_instantiation_of<SubExpression, simplify_t<E2>>>,
                std::negation<one_right<E1, E2>>,
                std::negation<zero_right<E1, E2>>
            >;

            template<typename E1, typename E2>
            using add_right = std::conjunction<
                internal::is_instantiation_of<AddExpression, simplify_t<E2>>,
                std::negation<one_left<E1, E2>>,
                std::negation<zero_left<E1, E2>>
            >;

            template<typename E1, typename E2>
            using add_left = std::conjunction<
                internal::is_instantiation_of<AddExpression, simplify_t<E1>>,
                std::negation<internal::is_instantiation_of<AddExpression, simplify_t<E2>>>,
                std::negation<sub_right<E1, E2>>,
                std::negation<one_right<E1, E2>>,
                std::negation<zero_right<E1, E2>>
            >;

            template<typename E1, typename E2>
            using same = std::conjunction<
                eq_t<simplify_t<E1>, simplify_t<E2>>,
                std::negation<internal::is_instantiation_of<ConstantExpression, simplify_t<E1>>>
            >;

            template<typename E1, typename E2>
            using sub_general = std::negation<
                std::disjunction<
                    two_constants<E1, E2>,
                    same<E1, E2>,
                    zero_right<E1, E2>,
                    zero_left<E1, E2>
                >
            >;

            template <typename E1, typename E2>
            using add_general = 
                std::negation<
                    std::disjunction<
                        zero_left<E1, E2>,
                        zero_right<E1, E2>,
                        two_constants<E1, E2>,
                        same<E1, E2>,
                        right_min_left<E1, E2>,
                        left_min_right<E1, E2>
                    >
                >;

            template <typename E1, typename E2>
            using mul_general = 
                std::negation<std::disjunction<
                    zero_left<E1, E2>,
                    zero_right<E1, E2>,
                    two_constants<E1, E2>,
                    one_left<E1, E2>,
                    one_right<E1, E2>,
                    sub_left<E1, E2>,
                    sub_right<E1, E2>,
                    add_right<E1, E2>,
                    add_left<E1, E2>,
                    two_exp<E1, E2>
                >>;
        };

        // e^0
        template<typename E>
        struct base_traits<ExpExpression<E>, 
                std::enable_if_t<eq_t<ZERO, simplify_t<E>>::value>
        > : base_traits<ONE> {};

        // x + y
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>, 
            std::enable_if_t<simplification_rules::template add_general<E1, E2>::value>>
        {
            using type = AddExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template add_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = AddExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = AddExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        // x - y
        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>, std::enable_if_t<
                    simplification_rules::template sub_general<E1, E2>::value
        >> {
            using type = SubExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template sub_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = SubExpression<typename base_traits<E1>::derive_t, typename base_traits<E2>::derive_t>;
            using derive_x = SubExpression<typename base_traits<E1>::derive_x, typename base_traits<E2>::derive_x>;
        };

        // -x
        template<typename E>
        struct base_traits<MinExpression<E>> {
            using type = MinExpression<E>;
            using value_at_t0 = FPQ64X::sub_t<typename FPQ64X::zero, typename base_traits<E>::value_at_t0>;
            using derive_t = MinExpression<typename base_traits<E>::derive_t>;
            using derive_x = MinExpression<typename base_traits<E>::derive_x>;
        };

        // x * y
        template<typename E1, typename E2>
        struct base_traits<
                MulExpression<E1, E2>, std::enable_if_t<simplification_rules::template mul_general<E1, E2>::value>
            > {
            using type = MulExpression<E1, E2>;
            using value_at_t0 = FPQ64X::template mul_t<typename base_traits<E1>::value_at_t0, typename base_traits<E2>::value_at_t0>;
            using derive_t = AddExpression<
                            MulExpression<typename base_traits<E1>::derive_t, E2>,
                            MulExpression<E1, typename base_traits<E2>::derive_t>>;
            
            using derive_x = AddExpression<
                            MulExpression<typename base_traits<E1>::derive_x, E2>,
                            MulExpression<E1, typename base_traits<E2>::derive_x>>;
        };
        
        // simplification rules
        // a + b
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template two_constants<E1, E2>::value>
        > : base_traits<ConstantExpression<Q64::add_t<typename simplify_t<E1>::v, typename simplify_t<E2>::v>>> {};

        // a - b
        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template two_constants<E1, E2>::value>
        > : base_traits<ConstantExpression<Q64::sub_t<typename simplify_t<E1>::v, typename simplify_t<E2>::v>>> {};

        // a * b
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template two_constants<E1, E2>::value>
        > : base_traits<ConstantExpression<Q64::mul_t<typename simplify_t<E1>::v, typename simplify_t<E2>::v>>> {};

        // e^x * e^y
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template two_exp<E1, E2>::value>
        > : base_traits<ExpExpression<AddExpression<typename simplify_t<E1>::v, typename simplify_t<E2>::v>>> {};

        // x + x
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template same<E1, E2>::value>
        > : base_traits<MulExpression<TWO, E2>> {};

        // x + (-x)
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                std::enable_if_t<simplification_rules::template right_min_left<E1, E2>::value>
        > : base_traits<ZERO> {};
        
        // -x + x
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                std::enable_if_t<simplification_rules::template left_min_right<E1, E2>::value>
        > : base_traits<ZERO> {};

        // x + 0 
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template zero_right<E1, E2>::value>
        > : base_traits<E1> {};

        // x * 0
         template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template zero_right<E1, E2>::value>
        > : base_traits<ZERO> {};

        // x - x
        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>, std::enable_if_t<simplification_rules::template same<E1, E2>::value>> : base_traits<ZERO> {};

        // x - 0
        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>, std::enable_if_t<simplification_rules::template zero_right<E1, E2>::value>> : base_traits<E1> {};

         // 0 - x
        template<typename E1, typename E2>
        struct base_traits<SubExpression<E1, E2>, std::enable_if_t<simplification_rules::template zero_left<E1, E2>::value>> : base_traits<MinExpression<E2>> {};

        // 0 + x
        template<typename E1, typename E2>
        struct base_traits<AddExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template zero_left<E1, E2>::value>
        > : base_traits<E2> {};
        
        // 0 * x
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template zero_left<E1, E2>::value>
        > : base_traits<ZERO> {};

        // 1 * x
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template one_left<E1, E2>::value>
        > : base_traits<E2> {};

        // x * 1
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>,
                    std::enable_if_t<simplification_rules::template one_right<E1, E2>::value>
        > : base_traits<E1> {};

        // x * (y - z)
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>, std::enable_if_t<simplification_rules::template sub_right<E1, E2>::value>
        > : base_traits<SubExpression<MulExpression<simplify_t<E1>, typename simplify_t<E2>::left>, MulExpression<simplify_t<E1>, typename simplify_t<E2>::right>>> {};

        // (x - y) * z
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>, std::enable_if_t<simplification_rules::template sub_left<E1, E2>::value>
        > : base_traits<SubExpression<MulExpression<typename simplify_t<E1>::left, simplify_t<E2>>, MulExpression<typename simplify_t<E1>::right, simplify_t<E2>>>> {};

        // x * (y + z)
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>, std::enable_if_t<simplification_rules::template add_right<E1, E2>::value>
        > : base_traits<AddExpression<MulExpression<simplify_t<E1>, typename simplify_t<E2>::left>, MulExpression<simplify_t<E1>, typename simplify_t<E2>::right>>> {};
        
        // (x + y) * z
        template<typename E1, typename E2>
        struct base_traits<MulExpression<E1, E2>, std::enable_if_t<simplification_rules::template add_left<E1, E2>::value>
        > : base_traits<AddExpression<MulExpression<typename simplify_t<E1>::left, simplify_t<E2>>, MulExpression<typename simplify_t<E1>::right, simplify_t<E2>>>> {};

        template<typename E, size_t k>
        struct taylor_coeff_helper {
            using type = FPQ64X::remove_divisor_t<FPQ64X::template div_t<typename taylor_coeff_helper<typename simplify_t<E>::derive_t, k-1>::type,
                                                FPQ64X::inject_constant_t<k>>>;
        };

        template<typename E>
        struct taylor_coeff_helper<E, 0> {
            using type = typename simplify_t<E>::value_at_t0;
        };

        template<typename E, size_t k>
        using taylor_coeff_t =  typename taylor_coeff_helper<E, k>::type::x;


        template<typename E, typename I>
        struct taylor_low {};

        template<typename E, size_t... k>
        struct taylor_low<E, std::index_sequence<k...>> {
            using type = P64T::template val<taylor_coeff_t<E, k>...>;
        };

        template<typename E, size_t k>
        using taylor_expansion_t = typename taylor_low<E, internal::make_index_sequence_reverse<k>>::type;
    }
}