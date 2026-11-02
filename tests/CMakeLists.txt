include(GoogleTest)

add_executable(lib_tests
  concepts_tests.cpp
  embedding_tests.cpp
  evaluation_tests.cpp
  fraction_field_tests.cpp
  integers_tests.cpp
  known_polynomials_tests.cpp
  polynomial_tests.cpp
  quotient_tests.cpp
  taylor_tests.cpp
  type_list_tests.cpp
  utilities_tests.cpp
  zpz_tests.cpp
)
target_link_libraries(lib_tests PRIVATE aerobus aerobus_dev_flags GTest::gtest_main)
gtest_discover_tests(lib_tests)

# install this build tree, then build a consumer through find_package(Aerobus)
add_test(NAME package.find_package
  COMMAND "${CMAKE_COMMAND}"
    "-DAEROBUS_BINARY_DIR=${PROJECT_BINARY_DIR}"
    "-DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/package_test"
    "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
    "-DCONFIG=$<CONFIG>"
    -P "${CMAKE_CURRENT_SOURCE_DIR}/package/check_package.cmake")
