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
instead of runtime. Using these capabilities, vastly extended with latest versions of the standard, we provide a library solving those issues. 
Polynomials can have their coefficients in any discrete integral domain. We provide a way to generate the fraction field of such rings -- including polynomials. 
And a way to generate taylor series of any math functions as long as coefficients are known. 

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

# Acknowledgements

Many thanks to my math teachers, A. Soyeur and M. Gonnord. I also acknowledge indirect contributions from F. Duguet, 
who basically learnt me all I know in C++. 

# Reference
