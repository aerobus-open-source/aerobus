#!/usr/bin/env python3
"""Differential test of aerobus::bigint against Python's integers.

Draws random operands (seeded, so that a failure reproduces), mixing random digits with digits that stress the
carries and the quotient estimate of Knuth's algorithm D (0, 1, 0x7fffffff, 0x80000000, 0xffffffff...), writes one
static_assert per result, and compiles them as one translation unit.

Usage: tools/check_bigint.py <c++ compiler> [--seed N] [--count N] [extra flags...]
"""
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASE = 1 << 32
SPECIAL = [0, 1, 2, 0x7FFFFFFF, 0x80000000, 0x80000001, 0xFFFFFFFE, 0xFFFFFFFF]


def as_type(x):
    negative = x < 0
    magnitude = -x if negative else x
    digits = []
    while magnitude:
        digits.append(magnitude & (BASE - 1))
        magnitude >>= 32
    return "V<" + ("true" if negative else "false") + "".join(f", 0x{d:x}u" for d in reversed(digits)) + ">"


def operand(rng, digits):
    x = sum((rng.choice(SPECIAL) if rng.random() < 0.5 else rng.getrandbits(32)) << (32 * i) for i in range(digits))
    return -x if rng.random() < 0.5 else x


def truncated_division(a, b):
    q = abs(a) // abs(b)
    q = q if (a < 0) == (b < 0) else -q
    return q, a - q * b


def program(seed, count):
    rng = random.Random(seed)
    lines = ['#include "aerobus.h"',
             "template<bool negative, uint32_t... digits>",
             "using V = aerobus::bigint::val<negative, digits...>;",
             "using aerobus::bigint;"]
    for _ in range(count):
        a, b = operand(rng, rng.randint(1, 12)), operand(rng, rng.randint(1, 8))
        if b != 0:
            q, r = truncated_division(a, b)
            lines.append(f"static_assert(std::is_same_v<bigint::div_t<{as_type(a)}, {as_type(b)}>, {as_type(q)}> && "
                         f"std::is_same_v<bigint::mod_t<{as_type(a)}, {as_type(b)}>, {as_type(r)}>);")
        a, b = operand(rng, rng.randint(0, 8)), operand(rng, rng.randint(0, 8))
        lines.append(f"static_assert(std::is_same_v<bigint::mul_t<{as_type(a)}, {as_type(b)}>, {as_type(a * b)}> && "
                     f"std::is_same_v<bigint::add_t<{as_type(a)}, {as_type(b)}>, {as_type(a + b)}> && "
                     f"std::is_same_v<bigint::sub_t<{as_type(a)}, {as_type(b)}>, {as_type(a - b)}>);")
    lines.append("int main() { return 0; }")
    return "\n".join(lines) + "\n", len(lines) - 5


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__, file=sys.stderr)
        return 2
    compiler, seed, count, flags = args[0], 2026, 1000, []
    rest = iter(args[1:])
    for arg in rest:
        if arg == "--seed":
            seed = int(next(rest))
        elif arg == "--count":
            count = int(next(rest))
        else:
            flags.append(arg)
    source_text, checks = program(seed, count)
    with tempfile.TemporaryDirectory() as work:
        source = Path(work, "bigint_differential.cpp")
        source.write_text(source_text)
        result = subprocess.run([compiler, "-std=c++20", "-fsyntax-only", f"-I{ROOT / 'src'}", *flags, str(source)],
                                capture_output=True, text=True)
    if result.returncode != 0:
        first = next((line for line in result.stderr.splitlines() if "error" in line), result.stderr)
        print(f"bigint differs from Python (seed {seed}, count {count}):\n  {first}")
        return 1
    print(f"{checks} static_asserts (seed {seed}): bigint agrees with Python")
    return 0


if __name__ == "__main__":
    sys.exit(main())
