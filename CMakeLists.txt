cmake_minimum_required(VERSION 3.14...3.31)

project(
    Aerobus
    VERSION 1.2
    DESCRIPTION "Header-only C++20 library for compile-time algebra on polynomials, discrete rings and fields"
    LANGUAGES CXX
)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(AEROBUS_IS_TOP_LEVEL ON)
else()
  set(AEROBUS_IS_TOP_LEVEL OFF)
endif()

option(AEROBUS_BUILD_TESTS "Build the unit tests" ${AEROBUS_IS_TOP_LEVEL})
option(AEROBUS_BUILD_BENCHMARKS "Build the benchmarks (needs Google Benchmark; uses OpenMP when available)" OFF)

if(AEROBUS_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# the library: a single header
include(GNUInstallDirs)
add_library(aerobus INTERFACE)
add_library(aerobus::aerobus ALIAS aerobus)
target_include_directories(aerobus INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(aerobus INTERFACE cxx_std_20)

# installation and find_package(Aerobus) support
include(CMakePackageConfigHelpers)
set(AEROBUS_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/cmake/Aerobus")
install(FILES src/aerobus.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS aerobus EXPORT AerobusTargets)
install(EXPORT AerobusTargets NAMESPACE aerobus:: DESTINATION "${AEROBUS_CMAKE_INSTALL_DIR}")
configure_package_config_file(cmake/AerobusConfig.cmake.in
  "${PROJECT_BINARY_DIR}/AerobusConfig.cmake"
  INSTALL_DESTINATION "${AEROBUS_CMAKE_INSTALL_DIR}")
write_basic_package_version_file("${PROJECT_BINARY_DIR}/AerobusConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion
  ARCH_INDEPENDENT)
install(FILES
  "${PROJECT_BINARY_DIR}/AerobusConfig.cmake"
  "${PROJECT_BINARY_DIR}/AerobusConfigVersion.cmake"
  DESTINATION "${AEROBUS_CMAKE_INSTALL_DIR}")

# warnings and limits shared by every target built in this repository
add_library(aerobus_dev_flags INTERFACE)
if(MSVC)
  target_compile_options(aerobus_dev_flags INTERFACE /W4 /permissive- /bigobj /we4309 /we4756)
  target_compile_definitions(aerobus_dev_flags INTERFACE _USE_MATH_DEFINES)
else()
  # deep compile-time recursions (Conway polynomials, high degree Taylor series) need a large template depth
  target_compile_options(aerobus_dev_flags INTERFACE
    -Wall -Werror=overflow -ftemplate-depth=30000 -ftemplate-backtrace-limit=0)
endif()

if(AEROBUS_BUILD_TESTS)
  include(FetchContent)
  if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
  endif()
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
    URL_HASH SHA256=40d4ec942217dcc84a9ebe2a68584ada7d4a33a8ee958755763278ea1c5e18ff
  )
  # use the same C runtime as the tests on Windows
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)

  enable_testing()
  include(GoogleTest)

  add_executable(lib_tests src/lib_tests.cpp)
  target_link_libraries(lib_tests PRIVATE aerobus aerobus_dev_flags GTest::gtest_main)
  gtest_discover_tests(lib_tests)

  # install this build tree, then build a consumer through find_package(Aerobus)
  add_test(NAME package.find_package
    COMMAND "${CMAKE_COMMAND}"
      "-DAEROBUS_BINARY_DIR=${PROJECT_BINARY_DIR}"
      "-DWORK_DIR=${PROJECT_BINARY_DIR}/package_test"
      "-DCXX_COMPILER=${CMAKE_CXX_COMPILER}"
      "-DCONFIG=$<CONFIG>"
      -P "${PROJECT_SOURCE_DIR}/tests/package/check_package.cmake")
endif()

if(AEROBUS_BUILD_BENCHMARKS)
  find_package(benchmark CONFIG QUIET)
  if(NOT benchmark_FOUND)
    include(FetchContent)
    FetchContent_Declare(
      googlebenchmark
      GIT_REPOSITORY https://github.com/google/benchmark.git
      GIT_TAG v1.9.1
    )
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
    set(BENCHMARK_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googlebenchmark)
  endif()
  find_package(OpenMP)

  set(AEROBUS_BENCHMARK_FLAGS "$<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-march=native>"
      CACHE STRING "Extra compile options for the benchmarks (instruction set)")

  add_executable(benchmarks src/benchmarks.cpp)
  target_link_libraries(benchmarks PRIVATE aerobus aerobus_dev_flags benchmark::benchmark)
  target_compile_options(benchmarks PRIVATE ${AEROBUS_BENCHMARK_FLAGS})
  if(OpenMP_CXX_FOUND)
    target_link_libraries(benchmarks PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
