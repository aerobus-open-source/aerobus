#!/usr/bin/env python3
"""Generate the Conway polynomial table of src/aerobus.h from imports/conwaypolynomials.txt.

The data is Frank Luebeck's list of Conway polynomials: [p, n, [c0, c1, ..., cn]], constant term first.
aerobus keeps the primes p < 1000 and, to bound compile time, the degrees n <= max_degree(p).

Usage:
    tools/generate_conway.py          rewrite the table in src/aerobus.h
    tools/generate_conway.py --check  exit 1 if src/aerobus.h differs from the generated table
"""
import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "imports" / "conwaypolynomials.txt"
HEADER = ROOT / "src" / "aerobus.h"
BEGIN = "    // NOLINTBEGIN generated table, see tools/\n"
END = "    // NOLINTEND\n"


def max_degree(p):
    """Largest degree n kept for the prime p."""
    if p < 23:
        return 20
    if p < 263:
        return 19
    if p < 307:
        return 12
    return 9


def table(data):
    lines = []
    for p, n, coeffs in re.findall(r"\[(\d+),(\d+),\[([0-9,]*)\]\]", data):
        p, n = int(p), int(n)
        if p >= 1000 or n > max_degree(p):
            continue
        values = ", ".join(f"ZPZV<{c}>" for c in reversed(coeffs.split(",")))
        lines.append(f"    template<> struct ConwayPolynomial<{p}, {n}> {{ using ZPZ = aerobus::zpz<{p}>; "
                     f"using type = POLYV<{values}>; }};\n")
    return "".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="only check that the header is up to date")
    args = parser.parse_args()

    header = HEADER.read_text()
    polyv = "    #define POLYV aerobus::polynomial<ZPZ>::template val\n"
    start = header.index(polyv, header.index(BEGIN)) + len(polyv)
    # the table ends at the #undef of the helper macros
    end = header.index("    #undef ZPZV\n", start)
    generated = table(DATA.read_text())
    if args.check:
        if header[start:end] != generated:
            print("src/aerobus.h: the Conway table differs from tools/generate_conway.py output", file=sys.stderr)
            return 1
        print("Conway table up to date")
        return 0
    HEADER.write_text(header[:start] + generated + header[end:])
    return 0


if __name__ == "__main__":
    sys.exit(main())
