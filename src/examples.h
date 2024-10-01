#ifndef SRC_EXAMPLES_H_
#define SRC_EXAMPLES_H_
/**
 * \example examples/hermite.cpp
 * How to use aerobus::known_polynomials::hermite_phys polynomials
 */


/**
 * \example examples/custom_taylor.cpp
 * How to implement your own Taylor serie using aerobus::taylor
 */

/**
 * \example examples/fp16.cu
 * How to leverage CUDA __half and __half2 16 bits floating points number using aerobus::i16
 * Warning : due to an NVIDIA bug (lack of constexpr operators), performance is not good
 */

/**
 * \example examples/continued_fractions.cpp
 * How to use aerobus::ContinuedFraction to get approximations of known numbers
 */

/**
 * \example examples/modular_arithmetic.cpp
 * How to use aerobus::zpz to perform computations on rational fractions with coefficients in modular rings
 */


/**
 * \example examples/make_polynomial.cpp
 * How to build your own sequence of known polynomials, here [Abel polynomials](https://en.wikipedia.org/wiki/Abel_polynomials)
 * 
 */

/**
 * \example examples/polynomials_over_finite_field.cpp
 * How to build a known polynomial (here aerobus::known_polynomials::allone) with coefficients in a finite field
 * (here aerobus::zpz<2>) and get its value when evaluated at a value in this field (here 1). 
 */
#endif  // SRC_EXAMPLES_H_
