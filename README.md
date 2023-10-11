# aerobus
Aerobus is a C++20 template library which expresses algebra concepts in types

It defines some concepts, such a Ring, IntegralDomain or Field which can be used to construct the field of fractions of a Ring, polynomials with coefficients in such a set or rational fractions, all at compile time. 

In addition, it provides predefined polynomials, such as taylor series of usual functions. Given polynomials can then be evaluated at floating point values, aerobus therefore define compile time versions of usual functions, and very efficient implementations for runtime. 

## examples
### pure compile time
Lets consider the folliwing program, featuring function exp - 1, with 13 64 bits coefficients
```c++
int main() {
	using V = aerobus::expm1<aerobus::i64, 13>;
	static constexpr double xx = V::eval(0.1);
    printf("%lf\n", xx);
}
```
V AND xx are completely compiled at compile time, yielding the following assembly (clang 17)

```assembly
.LCPI0_0:
        .quad   0x3fbaec7b35a00d3a              # double 0.10517091807564763
main:                                   # @main
        push    rax
        lea     rdi, [rip + .L.str]
        movsd   xmm0, qword ptr [rip + .LCPI0_0] # xmm0 = mem[0],zero
        mov     al, 1
        call    printf@PLT
        xor     eax, eax
        pop     rcx
        ret
.L.str:
        .asciz  "%lf\n"
```
### evaluations on variables
On the other hand, one might want to define a runtime function this way : 

```c++
double expm1(const double x) {
	using V = aerobus::expm1<aerobus::i64, 13>;
	return V::eval(x);
}
```
again, coefficients are all computed compile time, yielding following assembly (given processor supports fused multiply add) : 
```assembly
.LCPI0_0:
        .quad   0x3de6124613a86d09              # double 1.6059043836821613E-10
.LCPI0_1:
        .quad   0x3e21eed8eff8d898              # double 2.08767569878681E-9
.LCPI0_2:
        .quad   0x3e5ae64567f544e4              # double 2.505210838544172E-8
.LCPI0_3:
        .quad   0x3e927e4fb7789f5c              # double 2.7557319223985888E-7
.LCPI0_4:
        .quad   0x3ec71de3a556c734              # double 2.7557319223985893E-6
.LCPI0_5:
        .quad   0x3efa01a01a01a01a              # double 2.4801587301587302E-5
.LCPI0_6:
        .quad   0x3f2a01a01a01a01a              # double 1.9841269841269841E-4
.LCPI0_7:
        .quad   0x3f56c16c16c16c17              # double 0.0013888888888888889
.LCPI0_8:
        .quad   0x3f81111111111111              # double 0.0083333333333333332
.LCPI0_9:
        .quad   0x3fa5555555555555              # double 0.041666666666666664
.LCPI0_10:
        .quad   0x3fc5555555555555              # double 0.16666666666666666
.LCPI0_11:
        .quad   0x3fe0000000000000              # double 0.5
.LCPI0_12:
        .quad   0x3ff0000000000000              # double 1
expm1(double):                              # @expm1(double)
        vxorpd  xmm1, xmm1, xmm1
        vmovsd  xmm2, qword ptr [rip + .LCPI0_0] # xmm2 = mem[0],zero
        vfmadd231sd     xmm2, xmm0, xmm1        # xmm2 = (xmm0 * xmm1) + xmm2
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_1] # xmm2 = (xmm0 * xmm2) +   
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_2] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_3] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_4] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_5] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_6] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_7] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_8] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_9] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_10] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_11] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm2, xmm0, qword ptr [rip + .LCPI0_12] # xmm2 = (xmm0 * xmm2) + mem
        vfmadd213sd     xmm0, xmm2, xmm1        # xmm0 = (xmm2 * xmm0) + xmm1
        ret
```
### apply on vectors and get proper vectorization
If applied to a vector of data, with proper compiler hints, gcc can easily generate vectorized version of the code : 

```c++
double compute_expm1(const size_t N, const double* const __restrict in, double* const __restrict out) {
	using V = aerobus::expm1<aerobus::i64, 13>;
    for (size_t i = 0; i < N; ++i) {
        out[i] = V::eval(in[i]);
    }
}
```

yielding : 

```assembly
compute_expm1(unsigned long, double const*, double*):
        lea     rax, [rdi-1]
        cmp     rax, 2
        jbe     .L5
        mov     rcx, rdi
        xor     eax, eax
        vxorpd  xmm1, xmm1, xmm1
        vbroadcastsd    ymm14, QWORD PTR .LC1[rip]
        vbroadcastsd    ymm13, QWORD PTR .LC3[rip]
        shr     rcx, 2
        vbroadcastsd    ymm12, QWORD PTR .LC5[rip]
        vbroadcastsd    ymm11, QWORD PTR .LC7[rip]
        sal     rcx, 5
        vbroadcastsd    ymm10, QWORD PTR .LC9[rip]
        vbroadcastsd    ymm9, QWORD PTR .LC11[rip]
        vbroadcastsd    ymm8, QWORD PTR .LC13[rip]
        vbroadcastsd    ymm7, QWORD PTR .LC15[rip]
        vbroadcastsd    ymm6, QWORD PTR .LC17[rip]
        vbroadcastsd    ymm5, QWORD PTR .LC19[rip]
        vbroadcastsd    ymm4, QWORD PTR .LC21[rip]
        vbroadcastsd    ymm3, QWORD PTR .LC23[rip]
        vbroadcastsd    ymm2, QWORD PTR .LC25[rip]
.L3:
        vmovupd ymm15, YMMWORD PTR [rsi+rax]
        vmovapd ymm0, ymm15
        vfmadd132pd     ymm0, ymm14, ymm1
        vfmadd132pd     ymm0, ymm13, ymm15
        vfmadd132pd     ymm0, ymm12, ymm15
        vfmadd132pd     ymm0, ymm11, ymm15
        vfmadd132pd     ymm0, ymm10, ymm15
        vfmadd132pd     ymm0, ymm9, ymm15
        vfmadd132pd     ymm0, ymm8, ymm15
        vfmadd132pd     ymm0, ymm7, ymm15
        vfmadd132pd     ymm0, ymm6, ymm15
        vfmadd132pd     ymm0, ymm5, ymm15
        vfmadd132pd     ymm0, ymm4, ymm15
        vfmadd132pd     ymm0, ymm3, ymm15
        vfmadd132pd     ymm0, ymm2, ymm15
        vfmadd132pd     ymm0, ymm1, ymm15
        vmovupd YMMWORD PTR [rdx+rax], ymm0
        add     rax, 32
        cmp     rcx, rax
        jne     .L3
        mov     rax, rdi
        and     rax, -4
        vzeroupper
```

## predefined functions
```C++
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
```

## HOW TO
Clone or download the repo somewhere

include "lib.h" in code

Compile with -std=c++20 (at least)


[![DOI](https://zenodo.org/badge/499577459.svg)](https://zenodo.org/badge/latestdoi/499577459)
