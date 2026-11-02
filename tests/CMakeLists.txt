include(GoogleTest)

add_executable(lib_tests
  api_tests.cpp
  concepts_tests.cpp
  embedding_tests.cpp
  evaluation_tests.cpp
  fraction_field_tests.cpp
  integers_tests.cpp
  known_polynomials_tests.cpp
  macros_tests.cpp
  memory_tests.cpp
  polynomial_tests.cpp
  quotient_tests.cpp
  ring_axioms_tests.cpp
  taylor_tests.cpp
  type_list_tests.cpp
  utilities_tests.cpp
  zpz_tests.cpp
)
target_link_libraries(lib_tests PRIVATE aerobus aerobus_dev_flags GTest::gtest_main)
target_compile_definitions(lib_tests PRIVATE AEROBUS_CMAKE_PROJECT_VERSION="${PROJECT_VERSION}")
gtest_discover_tests(lib_tests)

# compile_fail/<name>.cpp must not compile, and the compiler output must match <regex>
# (the message of the static_assert that rejects it)
function(aerobus_compile_fail_test name regex)
  add_executable(compile_fail_${name} EXCLUDE_FROM_ALL compile_fail/${name}.cpp)
  target_link_libraries(compile_fail_${name} PRIVATE aerobus aerobus_dev_flags)
  add_test(NAME compile_fail.${name}
    COMMAND "${CMAKE_COMMAND}" --build "${PROJECT_BINARY_DIR}" --target compile_fail_${name} --config $<CONFIG>)
  set_tests_properties(compile_fail.${name} PROPERTIES
    PASS_REGULAR_EXPRESSION "${regex}"
    RESOURCE_LOCK aerobus_build_tree)
endfunction()

aerobus_compile_fail_test(zpz_division_by_zero "division by zero in zpz")
aerobus_compile_fail_test(zpz_non_invertible_divisor "divisor is not invertible modulo p")
aerobus_compile_fail_test(polynomial_division_not_exact "polynomial division is not exact")
aerobus_compile_fail_test(polynomial_division_by_zero "polynomial division by zero")
aerobus_compile_fail_test(polynomial_foreign_coefficient "coefficients of polynomial<Ring> must be values of Ring")
aerobus_compile_fail_test(fraction_zero_by_zero "division by zero in FractionField")
aerobus_compile_fail_test(fraction_one_by_zero "division by zero in FractionField")
aerobus_compile_fail_test(fraction_field_of_zero_divisors "FractionField needs an integral domain")
aerobus_compile_fail_test(conway_missing "no Conway polynomial for this")
aerobus_compile_fail_test(polynomial_horner_reduction_is_private "is private within this context|is a private member")

# install this build tree, then build a consumer through find_package(Aerobus)
add_test(NAME package.find_package
  COMMAND "${CMAKE_COMMAND}"
    "-DAEROBUS_BINARY_DIR=${PROJECT_BINARY_DIR}"
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/package_test"
    "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCONFIG=$<CONFIG>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/package/check_package.cmake")

# the header compiles without warnings when its templates are instantiated
if(NOT MSVC)
  add_executable(hygiene hygiene/instantiate.cpp)
  target_link_libraries(hygiene PRIVATE aerobus)
  target_compile_options(hygiene PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

# host-side checks of the CUDA fp16 code (WITH_CUDA_FP16), when nvcc is available; no GPU needed
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 75)
  endif()
  enable_language(CUDA)
  add_executable(cuda_fp16_host_tests cuda/fp16_host_tests.cu)
  target_link_libraries(cuda_fp16_host_tests PRIVATE aerobus)
  set_target_properties(cuda_fp16_host_tests PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON)
  target_compile_options(cuda_fp16_host_tests PRIVATE --expt-relaxed-constexpr)
  add_test(NAME cuda.fp16_host COMMAND cuda_fp16_host_tests)
endif()

# the Conway table of aerobus.h is the output of tools/generate_conway.py
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
  add_test(NAME conway.table_is_generated
    COMMAND "${Python3_EXECUTABLE}" -I "${PROJECT_SOURCE_DIR}/tools/generate_conway.py" --check)
endif()

# every C++ block of README.md compiles
if(Python3_Interpreter_FOUND AND NOT MSVC)
  add_test(NAME readme.cpp_blocks_compile
    COMMAND "${Python3_EXECUTABLE}" -I "${PROJECT_SOURCE_DIR}/tools/check_readme.py" "${CMAKE_CXX_COMPILER}")
endif()
