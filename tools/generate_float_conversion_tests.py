#!/usr/bin/env python3
"""Generate tests/float_conversion_tests.cpp: bigint and fractions converted to float, double and x87 long double
must be correctly rounded (to nearest, ties to even). The expected values come from exact fractions, with an
independent rounding checked against CPython's correctly rounded float(Fraction).
Also generate tests/cuda/half_conversion_tests.cu, the same for __half and __half2.

Usage: tools/generate_float_conversion_tests.py [--check]
  --check  compare with the files in the repository instead of writing them
"""
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / "tests" / "float_conversion_tests.cpp"
HALF_OUTPUT = OUTPUT.parent / "cuda" / "half_conversion_tests.cu"
HALF = (11, -14, 15)
FORMATS = {"float": (24, -126, 127), "double": (53, -1022, 1023), "long double": (64, -16382, 16383)}


def round_exact(value, precision, emin, emax):
    """value >= 0, rounded to nearest, ties to even: (m, k) with m 2^k the result, or None for an overflow."""
    if value == 0:
        return 0, 0
    exponent = value.numerator.bit_length() - value.denominator.bit_length()
    if Fraction(2) ** exponent > value:
        exponent -= 1
    k = max(exponent, emin) - precision + 1
    scaled = value / Fraction(2) ** k
    m = math.floor(scaled)
    rest = scaled - m
    if rest > Fraction(1, 2) or (rest == Fraction(1, 2) and m % 2 == 1):
        m += 1
    if m * Fraction(2) ** k >= Fraction(2) ** (emax + 1):
        return None
    return m, k


def literal(value, type_name):
    negative = value < 0
    precision, emin, emax = FORMATS[type_name]
    rounded = round_exact(abs(value), precision, emin, emax)
    sign = "-" if negative else ""
    if rounded is None:
        return f"{sign}std::numeric_limits<{type_name}>::infinity()"
    m, k = rounded
    suffix = {"float": "F", "double": "", "long double": "L"}[type_name]
    return f"{sign}0x{m:x}p{k:+d}{suffix}"


def bigint_type(value):
    magnitude, digits = abs(value), []
    while magnitude:
        digits.append(magnitude & 0xFFFFFFFF)
        magnitude >>= 32
    return "B<" + ("true" if value < 0 else "false") + "".join(f", 0x{d:x}U" for d in reversed(digits)) + ">"


def fraction_type(value):
    return f"QB::val<{bigint_type(value.numerator)}, {bigint_type(value.denominator)}>"


def cases():
    p2 = lambda e: Fraction(2) ** e  # noqa: E731
    out = [Fraction(1, math.factorial(n)) for n in range(2, 39)]
    out += [Fraction(2**100 + 2**47 + 1), Fraction(3**60), Fraction(7**40), Fraction(-(2**64 + 2**11 + 1))]
    for precision, emin, emax in (FORMATS["float"], FORMATS["double"]):
        one = p2(precision - 1)
        # exact midpoints in both directions, and just around them
        out += [one * 2 + 1, one * 2 + 3, -(one * 2 + 1), one * 2 + 1 + p2(-40), one * 2 + 1 - p2(-40)]
        low = p2(emin - precision + 1)  # smallest subnormal
        out += [low, low / 2, low * 3 / 2, low / 2 + low / 2**20, low / 2 - low / 2**20, -low, low / 3, low / 8]
        out += [p2(emin), p2(emin) - low, p2(emin) - low / 2]  # the frontier of the subnormals
        largest = (2 - p2(1 - precision)) * p2(emax)
        half_ulp = p2(emax - precision)
        out += [largest, largest + half_ulp, largest + half_ulp - low, p2(emax + 1), -p2(emax + 1)]
    rng = random.Random(2026)
    for _ in range(60):
        num = rng.getrandbits(rng.randint(1, 300)) + 1
        den = rng.getrandbits(rng.randint(1, 300)) + 1
        out.append(Fraction(num if rng.random() < 0.7 else -num, den))
    return out


def half_bits(value):
    """The binary16 encoding of value, correctly rounded: m 2^k is encoded as ((k + 24) << 10) + m."""
    if value == 0:
        return 0
    rounded = round_exact(abs(value), *HALF)
    magnitude = 0x7C00 if rounded is None else ((rounded[1] + 24) << 10) + rounded[0]
    return magnitude | (0x8000 if value < 0 else 0)


def round_half(value):
    """The value of half_bits(value)."""
    magnitude = half_bits(value) & 0x7FFF
    exponent = magnitude >> 10
    significand = magnitude & 0x3FF if exponent == 0 else (magnitude & 0x3FF) + 1024
    return (-1 if half_bits(value) & 0x8000 else 1) * significand * Fraction(2) ** (max(exponent, 1) - 25)


def half_cases():
    p2 = lambda e: Fraction(2) ** e  # noqa: E731
    precision, emin, emax = HALF
    one = p2(precision - 1)
    low = p2(emin - precision + 1)  # smallest subnormal
    largest = (2 - p2(1 - precision)) * p2(emax)
    half_ulp = p2(emax - precision)
    out = [Fraction(0), Fraction(1), Fraction(-1), Fraction(1, 2), Fraction(1, 3), Fraction(-1, 3), Fraction(355, 113)]
    out += [Fraction(1, math.factorial(n)) for n in range(2, 12)]
    # exact midpoints in both directions, just around them, and so close that rounding to float first would hit them
    out += [one * 2 + 1, one * 2 + 3, -(one * 2 + 1), one * 2 + 1 + p2(-40), one * 2 + 1 - p2(-40)]
    out += [1 + p2(-11), 1 + p2(-11) + p2(-40), 1 + p2(-11) - p2(-40), 1 + 3 * p2(-11) - p2(-40)]
    out += [low, low / 2, low * 3 / 2, low / 2 + low / 2**20, low / 2 - low / 2**20, -low, low / 3, low / 8]
    out += [p2(emin), p2(emin) - low, p2(emin) - low / 2]  # the frontier of the subnormals
    out += [largest, largest + half_ulp, largest + half_ulp - low, p2(emax + 1), -p2(emax + 1), p2(200)]
    rng = random.Random(2026)
    for _ in range(60):
        num = rng.getrandbits(rng.randint(1, 120)) + 1
        den = rng.getrandbits(rng.randint(1, 120)) + 1
        out.append(Fraction(num if rng.random() < 0.7 else -num, den))
    return out


def half_program():
    w = ["// Generated by tools/generate_float_conversion_tests.py from exact fractions; do not edit by hand.",
         "// Built by nvcc and hipcc, runs without a GPU. The conversions are evaluated at compile time and read at run time,",
         "// std::bit_cast being not constexpr on the __half of HIP; the kernels show that they compile for the device.",
         "#include <bit>\n#include <cstdint>\n#include <cstdio>\n\n#include \"aerobus.h\"\n",
         "namespace {\ntemplate<bool negative, uint32_t... digits>\nusing B = aerobus::bigint::val<negative, digits...>;",
         "using QB = aerobus::FractionField<aerobus::bigint>;\n\n"
         "int failures = 0;\n\n"
         "void check(bool ok, int line) {\n    if (!ok) {\n"
         "        std::printf(\"FAILED: line %d\\n\", line);\n        ++failures;\n    }\n}\n\n"
         "uint16_t bits(__half x) {\n    return std::bit_cast<uint16_t>(x);\n}\n\n"
         "template<typename V>\nvoid half_is(uint16_t expected, int line) {\n"
         "    static constexpr __half value = V::template get<__half>();\n"
         "    check(bits(value) == expected, line);\n}\n\n"
         "template<typename V>\nvoid half2_is(uint16_t expected, int line) {\n"
         "    const __half2 value = V::template get<__half2>();\n"
         "    check(bits(__low2half(value)) == expected && bits(__high2half(value)) == expected, line);\n}\n"
         "}  // namespace\n"]
    cases = half_cases()
    for v in cases:
        # decoding the encoding gives a value within half a unit in the last place
        if half_bits(v) & 0x7FFF < 0x7C00:
            assert abs(round_half(v) - v) <= Fraction(2) ** (max(half_bits(v) >> 10 & 0x1F, 1) - 26), v
    w.append("""// the coefficients of a polynomial of tools/polynomial.sollya, evaluated on the device
using R = aerobus::polynomial<QB>::val<
    aerobus::from_hex_float_t<"0x1.0bp-7">,
    aerobus::from_hex_float_t<"-0x1.554p-3">,
    aerobus::from_hex_float_t<"0x1p0">>;

__global__ void half_kernel(const __half* x, __half* y) {
    y[0] = R::eval(x[0]);
    y[1] = B<false, 3U>::get<__half>();
}

__global__ void half2_kernel(const __half2* x, __half2* y) {
    y[0] = R::eval(x[0]);
    y[1] = B<false, 3U>::get<__half2>();
}

int main() {""")
    w.append("    // bigint")
    for v in cases:
        if v.denominator == 1:
            w.append(f"    half_is<{bigint_type(v.numerator)}>(0x{half_bits(v):04x}, __LINE__);")
    w.append("\n    // fractions of bigint")
    for v in cases:
        w.append(f"    half_is<{fraction_type(v)}>(0x{half_bits(v):04x}, __LINE__);")
    w.append("\n    // fractions of native integers")
    native = [Fraction(1, 3), Fraction(-1, 3), Fraction(355, 113), Fraction(2**11 + 1), Fraction(2**24 + 1, 2**13),
              Fraction(2**24 + 2**13 + 1, 2**24), Fraction(2**24 + 3 * 2**13 - 1, 2**24), Fraction(-(2**25 + 3), 7)]
    for v in native:
        t = f"aerobus::make_q32_t<{v.numerator}, {v.denominator}>"
        w.append(f"    half_is<{t}>(0x{half_bits(v):04x}, __LINE__);")
    w.append("\n    // __half2")
    for v in [Fraction(1, 3), Fraction(-1, 3), Fraction(1, 2), 1 + Fraction(2) ** -11 + Fraction(2) ** -40]:
        w.append(f"    half2_is<{fraction_type(v)}>(0x{half_bits(v):04x}, __LINE__);")
    w.append(f"    half2_is<{bigint_type(3)}>(0x{half_bits(Fraction(3)):04x}, __LINE__);")
    # the host evaluation of the polynomial above: Horner, with one rounding per fused multiply-add
    x = Fraction(1, 2)
    p = round_half(round_half(Fraction(0x10b, 2**15) * x + Fraction(-0x1554, 2**15)) * x + 1)
    w.append(f"""
    check(bits(R::eval(__float2half(0.5F))) == 0x{half_bits(p):04x}, __LINE__);
    std::printf("%d failure(s)\\n", failures);
    return failures == 0 ? 0 : 1;
}}""")
    w = [line + "  // NOLINT(whitespace/line_length)" if len(line) > 120 and "\n" not in line else line for line in w]
    return "\n".join(w) + "\n"


def write_or_check(path, text):
    if "--check" in sys.argv[1:]:
        if path.read_text() != text:
            print(f"{path} differs from the output of {Path(__file__).name}: regenerate it")
            return 1
        print(f"{path.name} is up to date")
        return 0
    path.write_text(text)
    print(f"wrote {path} ({len(text.splitlines())} lines)")
    return 0


def main():
    w = []
    w.append("// Generated by tools/generate_float_conversion_tests.py from exact fractions; do not edit by hand.")
    w.append("#include <gtest/gtest.h>\n#include <bit>\n#include <cfloat>\n#include <cstdint>\n#include <limits>\n\n#include \"aerobus.h\"\n")
    w.append("namespace {\ntemplate<bool negative, uint32_t... digits>\nusing B = aerobus::bigint::val<negative, digits...>;")
    w.append("using QB = aerobus::FractionField<aerobus::bigint>;\n}  // namespace\n")
    all_cases = cases()
    for value in all_cases:
        if abs(value) > 0:
            # the independent rounding agrees with CPython's correctly rounded conversion
            sign = -1.0 if value < 0 else 1.0
            try:
                expected = float(value)
            except OverflowError:
                expected = sign * math.inf
            rounded = round_exact(abs(value), *FORMATS["double"])
            mine = math.inf if rounded is None else float(rounded[0] * Fraction(2) ** rounded[1])
            assert sign * mine == expected, value
    integers = [v for v in all_cases if v.denominator == 1]
    for type_name in ("float", "double"):
        name = type_name.replace(" ", "_")
        w.append(f"TEST(float_conversion, bigint_to_{name}) {{")
        for v in integers:
            w.append(f"    static_assert({bigint_type(v.numerator)}::get<{type_name}>() == {literal(v, type_name)});")
        w.append("}\n")
        w.append(f"TEST(float_conversion, bigint_fractions_to_{name}) {{")
        for v in all_cases:
            w.append(f"    static_assert({fraction_type(v)}::get<{type_name}>() == {literal(v, type_name)});")
        w.append("}\n")
    w.append("#if LDBL_MANT_DIG == 64\nTEST(float_conversion, x87_long_double) {")
    for v in all_cases[:37] + integers[:8]:
        w.append(f"    static_assert({fraction_type(v)}::get<long double>() == {literal(v, 'long double')});")
    w.append("}\n#endif\n")
    # fractions of native integers
    q64 = [Fraction(1, math.factorial(n)) for n in range(2, 21)] + [Fraction(2**53 + 1), Fraction(2**62 + 1, 3)]
    q64 += [Fraction(10**18 + 9, 10**18 + 3), Fraction(-(2**53 + 3), 2**61 - 1)]
    w.append("TEST(float_conversion, q64_to_double_and_float) {")
    for v in q64:
        t = f"aerobus::make_q64_t<{v.numerator}LL, {v.denominator}LL>"
        w.append(f"    static_assert({t}::get<double>() == {literal(v, 'double')});")
        w.append(f"    static_assert({t}::get<float>() == {literal(v, 'float')});")
    w.append("}\n")
    q32 = [Fraction(1, 3), Fraction(100000001, 100000007), Fraction(2**24 + 1), Fraction(-(2**25 + 3), 7)]
    w.append("TEST(float_conversion, q32_to_float) {")
    for v in q32:
        w.append(f"    static_assert(aerobus::make_q32_t<{v.numerator}, {v.denominator}>::get<float>() == {literal(v, 'float')});")
    w.append("}\n")
    w.append("""TEST(float_conversion, zero_and_sign_of_zero) {
    static_assert(aerobus::bigint::zero::get<double>() == 0.0);
    static_assert(QB::zero::get<float>() == 0.0F);
    // a negative value rounded to zero keeps its sign
    using tiny = QB::val<B<true, 1U>, B<false, 1U, 0U, 0U, 0U, 0U, 0U>>;  // -2^-160
    static_assert(std::bit_cast<uint32_t>(tiny::get<float>()) == 0x80000000U);
}

TEST(float_conversion, integer_targets_are_unchanged) {
    // a conversion to an integer type still truncates the quotient
    static_assert(aerobus::make_q32_t<7, 2>::get<int>() == 3);
    static_assert(aerobus::make_q64_t<-7, 2>::get<int64_t>() == -3);
    static_assert(B<true, 12U>::get<int64_t>() == -12);
}""")
    w = [line + "  // NOLINT(whitespace/line_length)" if len(line) > 120 and "\n" not in line else line for line in w]
    text = "\n".join(w) + "\n"
    return write_or_check(OUTPUT, text) | write_or_check(HALF_OUTPUT, half_program())


if __name__ == "__main__":
    sys.exit(main())
