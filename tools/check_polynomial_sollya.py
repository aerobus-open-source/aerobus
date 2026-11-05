#!/usr/bin/env python3
"""Check tools/polynomial.sollya end to end: run it, compile its output, and compare the compiled polynomial with
the coefficients it printed, read here as exact fractions. Exits with 77 (skipped) without Sollya.

Usage: tools/check_polynomial_sollya.py <c++ compiler> [extra flags...]
"""
import math
import re
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "tools" / "polynomial.sollya"
SKIPPED = 77
ZERO = "FractionField<aerobus::bigint>::zero"
COEFFICIENT = re.compile(rf'from_hex_float_t<"(-?)0x([0-9a-f]+)(?:\.([0-9a-f]+))?p([+-]?[0-9]+)">|({ZERO})')
PRECISIONS = {"HP": 11, "SG": 24, "D": 53, "QD": 113, "DD": None, "TD": None}


def sinc(x):
    return math.sin(x) / x if x else 1.0


# Sollya arguments (None for the defaults), target and interval as Python sees them, bound on the error, degree
CASES = [
    (None, sinc, (-math.pi / 4, math.pi / 4), 1e-6, 24),
    (("sin(x)/x", "[-pi/4, pi/4]", "[|0, 2, 4, 6, 8|]", "SG"), sinc, (-math.pi / 4, math.pi / 4), 1e-6, 8),
    (("exp(x)", "[0, 1]", "7", "D"), math.exp, (0, 1), 1e-6, 7),
    (("cos(x)", "[-1/2, 1/2]", "[|0, 2, 4, 6|]", "HP"), math.cos, (-0.5, 0.5), 1e-6, 6),
    (("log(1 + x)", "[0, 1/2]", "[|1, 2, 3, 4, 5, 6|]", "QD"), math.log1p, (0, 0.5), 1e-6, 6),
    (("sin(x)", "[-1/2, 1/2]", "[|1, 3, 5, 7|]", "TD"), math.sin, (-0.5, 0.5), 1e-10, 7),
    (("sin(x)/x", "[-pi/4, pi/4]", "[|0, 2, 4, 6, 8, 10|]", "DD"), sinc, (-math.pi / 4, math.pi / 4), 1e-14, 10),
    (("exp(x)", "[-1/4, 1/4]", "", "SG"), math.exp, (-0.25, 0.25), 1e-6, None),  # guessed from the format
]


def run_script(arguments):
    command = ["sollya", "--nocolor", "--warnonstderr", str(SCRIPT)]
    if arguments is not None:
        command += ["--args", *arguments]
    return subprocess.run(command, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=300,
                          check=True).stdout


def coefficient(match):
    """The value of one printed coefficient, computed without aerobus."""
    if match[5]:
        return Fraction(0)
    fraction = match[3] or ""
    value = Fraction(int(match[2] + fraction, 16)) * Fraction(2) ** (int(match[4]) - 4 * len(fraction))
    return -value if match[1] else value


def check_format(coefficients, precision):
    """The coefficients are numbers of the format, not approximations of them."""
    for c in coefficients:
        assert c.denominator & (c.denominator - 1) == 0, c
        magnitude = abs(c.numerator)
        if magnitude:
            assert (magnitude // (magnitude & -magnitude)).bit_length() <= precision, c


def check_error(coefficients, target, interval, bound):
    low, high = interval
    for k in range(2001):
        x = Fraction(low) + (Fraction(high) - Fraction(low)) * k / 2000
        value = Fraction(0)
        for c in coefficients:
            value = value * x + c
        error = abs(float(value) - target(float(x)))
        assert error <= bound, f"error {error} above {bound} at {float(x)}"


def bigint_type(value):
    magnitude, digits = abs(value), []
    while magnitude:
        digits.append(magnitude & 0xFFFFFFFF)
        magnitude >>= 32
    sign = "true" if value < 0 else "false"
    return f"aerobus::bigint::val<{sign}" + "".join(f", 0x{d:x}U" for d in reversed(digits)) + ">"


def fraction_type(value):
    return f"Q::val<{bigint_type(value.numerator)}, {bigint_type(value.denominator)}>"


def program(outputs):
    lines = ['#include "aerobus.h"', "#include <type_traits>", "using Q = aerobus::FractionField<aerobus::bigint>;"]
    for index, (text, coefficients) in enumerate(outputs):
        lines += [f"namespace case{index} {{", text]
        lines.append(f"static_assert(type::degree == {len(coefficients) - 1});")
        for degree, c in enumerate(reversed(coefficients)):
            lines.append(f"static_assert(std::is_same_v<type::coeff_at_t<{degree}>, {fraction_type(c)}>);")
        lines.append(f"static_assert(type::coeff_at_t<0>::get<double>() == {float(coefficients[-1]).hex()});")
        lines.append("}  // namespace")
    return "\n".join(lines) + "\n"


def check_failure(arguments, expected):
    """A failure prints #error, which stops the compilation of the output, and no type."""
    text = run_script(arguments)
    assert expected in text and "using type" not in text, f"{arguments}: {text}"


def main():
    if shutil.which("sollya") is None:
        print("sollya is not installed")
        return SKIPPED
    compiler, flags = sys.argv[1], sys.argv[2:]
    check_failure(("x", "[0, 1]", "3", "ZZ"), "#error FORMAT")
    check_failure(("cos(x)", "[-1/2, 1/2]", "[|0, 2, 4|]", "SG"), "#error fpminimax")
    check_failure(("x",), "#error expected")
    outputs = []
    for arguments, target, interval, bound, degree in CASES:
        text = run_script(arguments)
        coefficients = [coefficient(m) for m in COEFFICIENT.finditer(text)]
        assert coefficients and "#error" not in text, f"{arguments}: {text}"
        assert degree is None or len(coefficients) == degree + 1, f"{arguments}: {text}"
        precision = PRECISIONS[arguments[3] if arguments else "SG"]
        if precision is not None:
            check_format(coefficients, precision)
        check_error(coefficients, target, interval, bound)
        outputs.append((text, coefficients))
    with tempfile.TemporaryDirectory() as directory:
        source = Path(directory) / "polynomial_sollya.cpp"
        source.write_text(program(outputs))
        result = subprocess.run([compiler, "-std=c++20", "-fsyntax-only", f"-I{ROOT / 'src'}", *flags, str(source)],
                                capture_output=True, text=True, check=False)
        if result.returncode != 0:
            print(result.stderr[:4000])
            return 1
    print(f"{len(outputs)} polynomials compiled and agree with the coefficients printed by Sollya")
    return 0


if __name__ == "__main__":
    sys.exit(main())
