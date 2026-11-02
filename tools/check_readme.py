#!/usr/bin/env python3
"""Compile every ```cpp block of README.md against src/aerobus.h.

Each block is compiled on its own, inside a namespace, with AEROBUS_CONWAY_IMPORTS defined.

Usage: tools/check_readme.py <c++ compiler> [extra flags...]
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    if len(sys.argv) < 2:
        print(__doc__, file=sys.stderr)
        return 2
    compiler, flags = sys.argv[1], sys.argv[2:]
    blocks = re.findall(r"```cpp\n(.*?)```", (ROOT / "README.md").read_text(), flags=re.S)
    failures = 0
    with tempfile.TemporaryDirectory() as work:
        for i, block in enumerate(blocks):
            source = Path(work, f"snippet_{i}.cpp")
            source.write_text("#define AEROBUS_CONWAY_IMPORTS\n#include \"aerobus.h\"\n"
                              f"namespace snippet_{i} {{\n{block}}}\nint main() {{ return 0; }}\n")
            result = subprocess.run([compiler, "-std=c++20", "-fsyntax-only", "-Wall", "-Werror",
                                     "-Wno-unused-variable", "-Wno-unused-const-variable", "-ftemplate-depth=30000", f"-I{ROOT / 'src'}", *flags, str(source)],
                                    capture_output=True, text=True)
            if result.returncode != 0:
                failures += 1
                first = next((line for line in result.stderr.splitlines() if "error" in line), result.stderr)
                print(f"README.md cpp block {i + 1} does not compile:\n{block}\n  {first}\n")
    print(f"{len(blocks) - failures}/{len(blocks)} README cpp blocks compile")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
