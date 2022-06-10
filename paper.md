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
Templates are a language-in-the-language which is Turing-complete, meaning we can run every computation at compile time
instead of runtime, as long as input data is known at compile time. 
Using these capabilities, vastly extended with latest versions of the standard, we implemented a library 
for polynomials over any discrete integral domain, such $\mathbb{Z}$. We  also provide a way to generate the fraction field of such 
rings (e.g. $\mathbb{Q}$) -- including polynomials. In addition, we expose a way to generate taylor series of any math functions 
as long as coefficients are known. 

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

# Mathematic definitions
A `ring` $\mathbb{A}$ is a non empty set with two internal laws, addition and multiplication. There is a neutral element for both, zero and one. 
Addition is commutative and associative and every element $x$ has an inverse $-x$. Multiplication is commutative and associative
and distributive over addition, meaning that $a(b+c) = ab+ac$ for every $a, b, c$ element. We call it `discrete` if it is countable. 

An `integral domain` is a ring with one additional property. For every elements $a, b, c$ such as $ab = ac$, then either $a = 0$ 
or $b = c$. Such a ring is not always a field, such as $\mathbb{Z}$ shows it. 

For such an integral domain, we can build two important structures : 
- the polynomials $\mathbb{A}[X]$, wich is in turn an integral domain when given usual operations ;
- the fraction field, which is the smallest field containg $\mathbb{A}$. Canonical example is $\mathbb{Q}$, the set of 
rational numbers. 

# Software


# Acknowledgements

Many thanks to my math teachers, A. Soyeur and M. Gonnord. I also acknowledge indirect contributions from F. Duguet, 
who basically learnt me all I know in C++. 

# Reference
