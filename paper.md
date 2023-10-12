---
title: 'Aerobus: a C++ template library for polynomials algebra over discrete integral domains'
tags:
  - Polynomials
  - Mathematics
  - metaprogramming
  - Integral Domains
  - Rings
  - Fields
authors:
  - name: Regis Portalez
    orcid: 0000-0002-6291-3345 
    affiliation: 1
affiliations:
 - name: Independant Researcher, France
   index: 1
date: 10 June 2022
bibliography: paper.bib
---

# Summary

C++ comes with high compile-time computations capability, also known as metaprogramming with templates.
Templates are a language-in-the-language which is Turing-complete, meaning we can run every computation at compile time instead of runtime, as long as input data is known at compile time. 

Using these capabilities, vastly extended with latest versions of the standard, we implemented a library 
for polynomials over any discrete integral domain, such $\mathbb{Z}$. We  also provide a way to generate the fraction field of such rings (e.g. $\mathbb{Q}$). 

We also implemented polynomials over such a discrete ring or field (e.g. $\mathbb{Q}[X]$). Since polynomials are also a ring, same implementation as above gives us rational fractions as field of fractions of polynomials.

In addition, we expose a way to generate taylor series of any math functions as long as coefficients are known. 

`Aerobus` was designed to be used in high performance software, teaching purposes or embedded sotfware where as much as possible must be precomputed to shrink binary size. It compiles with major compilers : gcc, clang and msvc. It is quite easily configurable and extensible. 

# Statement of need
Polynomials over discrete integral domains are highly used in physics simulations or cryptographic applications. 
Most of software applications evaluate transcendental functions, such as `exp`, `sin` by using C math library
which lead to a call -- inlining being often prevented by compiler when high precision is required. 
Hardwarde vendors provide high level libraries such as [@wang2014intel], but implementation is often hidden. 
We can also provide inline functions, such as [@vml] does. 
But library such VML are highly tight to one architecture by their use of intrinsics or inline assembly. 
In addition, they only provide a restricted list of math functions and do not expose capabilities to generate 
high performance versions of other functions such as arctanh. 
`Aerobus` provides automatic generation of such functions.
In addition, `Aerobus` provides a way to control precision of generated function by changing the degree of taylor expansion, which can't be use in competing libraries without reimplementing the whole function. 

# Mathematic definitions
A `ring` $\mathbb{A}$ is a non empty set with two internal laws, addition and multiplication. There is a neutral element for both, zero and one. 
Addition is commutative and associative and every element $x$ has an inverse $-x$. Multiplication is commutative, associative and distributive over addition, meaning that $a(b+c) = ab+ac$ for every $a, b, c$ element. We call it `discrete` if it is countable. 

An `integral domain` is a ring with one additional property. For every elements $a, b, c$ such as $ab = ac$, then either $a = 0$ or $b = c$. Such a ring is not always a field, such as $\mathbb{Z}$ shows it. 

For such an integral domain, we can build two important structures : 

- the polynomials $\mathbb{A}[X]$, wich is in turn an integral domain when given usual operations ;
- the field of fractions, which is the smallest field containg $\mathbb{A}$. Canonical example is $\mathbb{Q}$, the set of rational numbers. 

# Software
Library exposes a some native types: 

- `i32` and `i64` ($\mathbb{Z}$ seen as 32bits or 64 bits integers)
- `zpz` the quotient ring $\mathbb{Z}/p\mathbb{Z}$ 
- `polynomial<T>` where T is a ring
- `FractionField<T>` where T is an integral domain

Polynomial expose an evaluation function, which automatically generate Horner development and unrolls the
loop by generating it at compile time.

All of those types expose internal operations, such as: 

- add_t (addition)
- sub_t (substraction)
- mul_t (multiplication)
- is_zero_t (is zero element)
- one (one element -- neutral for multiplication)
- zero (zero elemet -- neutral for addition)
- eq_t (equality)

Integral domains and polynomials also expose additional operations : 

- div_t (division)
- mod_t (modulus)
- gcd_t (greatest common divisor)
- lt_t (less than)
- gt_t (greater than)

Library provide builtin integers and functions, such as:

- is_prime
- factorial_t
- pow_t
- alternate_t ($(-1)^p$)
- combination_t
- bernouilli_t

And taylor series for these functions:

- exp
- expm1 (exp - 1)
- lnp1 (ln(x+1))
- geom (1/(1-x))
- sin
- cos
- tan
- sh
- cosh
- tanh
- asin
- acos
- acosh
- asinh
- atanh

# Example and resulting assembly

Considering this function, applying expm1 twelve time to a double precision float: 
```
double expm1_12(const double x) {
	using V = aerobus::expm1<aerobus::i64, 13>;
	return V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(
		V::eval(V::eval(V::eval(V::eval(V::eval(V::eval(x))))))))))));
}

void vexpm1_12(const std::vector<double>& in, std::vector<double>& out) {
	for (int i = 0; i < in.size(); ++i) {
		out[i] = expm1_12(in[i]);
	}
}
```

Compiled with `clang 17.0` with options `-O3 -std=c++20 -mavx512f`, it yields an assembly similar to: 

```
vfmadd213pd %zmm0, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm0
vfmadd213pd %zmm2, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm2
vfmadd213pd %zmm3, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm3
vfmadd213pd %zmm4, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm4
vfmadd213pd %zmm5, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm5
vfmadd213pd %zmm6, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm6
vfmadd213pd %zmm7, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm7
vfmadd213pd %zmm8, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm8
vfmadd213pd %zmm9, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm9
vfmadd213pd %zmm10, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm10
vfmadd213pd %zmm11, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm11
vfmadd213pd %zmm12, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm12
vfmadd213pd %zmm13, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm13
vfmadd213pd %zmm1, %zmm15, %zmm14 # zmm14 = (zmm15 * zmm14) + zmm1
vxorpd %xmm15, %xmm15, %xmm15
```
In which all vector registers are saturated with fused multiply-add instructions. 

# Acknowledgements

Many thanks to my math teachers, A. Soyeur and M. Gonnord. I also acknowledge indirect contributions from F. Duguet, who basically learnt me all I know in C++. 

# Reference