namespace aerobus {
    struct Tchebychev {
        private:
            using T = aad::TExpression<Q64::one, 1>;
            using X = aad::XExpression<Q64::one, 1>;
            using XT = aad::MulExpression<X, T>;
            using ONE = aad::ConstantExpression<Q64::one>;
            using TWO = aad::ConstantExpression<Q64::inject_constant_t<2>>;
            using B0 = aad::DivExpression<
                // 1 - xt
                aad::SubExpression<ONE, XT>,
                // 1 - 2xt + t^2
                aad::AddExpression<
                    // 1 - 2xt
                    aad::SubExpression<ONE, aad::MulExpression<TWO, XT>>,
                    aad::TExpression<Q64::one, 2>
                >
            >;

        public:
            template<size_t k>
            using all = aad::taylor_expansion_t<B0, k>;
            
            template<size_t k>
            using at = aad::taylor_coeff_t<B0, k>;
    };
}