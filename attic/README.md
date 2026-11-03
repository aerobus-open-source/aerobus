# attic

Files kept for reference only. They are not built, tested, linted or documented.

They target `aerobus::libm` (sin, cos and friends evaluated with Sollya-generated polynomials), which was removed from `src/aerobus.h` in commit 18b714f. None of them compiles against the current header.

| File | Former location | Content |
| --- | --- | --- |
| `libm_tests.cpp` | `src/libm_tests.cpp` | Range-reduction and accuracy tests for `aerobus::libm::sin` / `cos` |
| `libm_benchmarks.cpp` | `src/benchmarks.cpp` | Google Benchmark cases for `libm::sin`, `libm::fast_sin`, `libm::cos` |
| `libm.cu` | `examples/libm.cu` | CUDA `__half` sine kernel |
| `libm_hip.cpp` | `examples/libm_hip.cpp` | HIP sine kernel |
| `sollya/script.sollya` | `src/sollya/script.sollya` | Sollya script emitting aerobus rational polynomials |
