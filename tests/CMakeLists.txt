include(GoogleTest)

add_executable(lib_tests
  concepts_tests.cpp
  embedding_tests.cpp
  evaluation_tests.cpp
  fraction_field_tests.cpp
  integers_tests.cpp
  known_polynomials_tests.cpp
  memory_tests.cpp
  polynomial_tests.cpp
  quotient_tests.cpp
  taylor_tests.cpp
  type_list_tests.cpp
  utilities_tests.cpp
  zpz_tests.cpp
)
target_link_libraries(lib_tests PRIVATE aerobus aerobus_dev_flags GTest::gtest_main)
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

# install this build tree, then build a consumer through find_package(Aerobus)
add_test(NAME package.find_package
  COMMAND "${CMAKE_COMMAND}"
    "-DAEROBUS_BINARY_DIR=${PROJECT_BINARY_DIR}"
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/package_test"
    "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCONFIG=$<CONFIG>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/package/check_package.cmake")
