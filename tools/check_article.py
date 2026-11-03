#!/usr/bin/env python3
"""Compile the ```cpp blocks of an article as one translation unit against src/aerobus.h.

The blocks of an article build on each other, so they are concatenated in order, after the standard
headers they use and with AEROBUS_CONWAY_IMPORTS defined. When <article>.check.cpp exists, it is
appended: static_asserts there check the numbers quoted in the text.

Usage: tools/check_article.py <article.md> <c++ compiler> [extra flags...]
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PRELUDE = """#include <array>
#include <cstddef>
#include <cstdint>
#include <utility>
#define AEROBUS_CONWAY_IMPORTS
#include "aerobus.h"
using std::size_t;
"""


def main():
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    article, compiler, flags = Path(sys.argv[1]), sys.argv[2], sys.argv[3:]
    blocks = re.findall(r"```cpp\n(.*?)```", article.read_text(), flags=re.S)
    if not blocks:
        print(f"{article}: no cpp block")
        return 1
    with tempfile.TemporaryDirectory() as work:
        source = Path(work, "article.cpp")
        # optional checks of the numbers quoted in the article, next to it as <article>.check.cpp
        checks = article.with_suffix(".check.cpp")
        epilogue = checks.read_text() if checks.exists() else ""
        source.write_text(PRELUDE + "\n".join(blocks) + "\n" + epilogue + "\nint main() { return 0; }\n")
        result = subprocess.run([compiler, "-std=c++20", "-fsyntax-only", "-Wall", "-Werror", "-Wno-unused-variable",
                                 "-Wno-unused-const-variable", "-ftemplate-depth=30000", f"-I{ROOT / 'src'}", *flags,
                                 str(source)], capture_output=True, text=True)
    if result.returncode != 0:
        print(f"{article}: the cpp blocks do not compile\n{result.stderr}")
        return 1
    print(f"{article}: {len(blocks)} cpp blocks compile")
    return 0


if __name__ == "__main__":
    sys.exit(main())
