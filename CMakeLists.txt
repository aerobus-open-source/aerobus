cmake_minimum_required(VERSION 3.1...3.27)

project(
    Aerobus
    VERSION 1.2
    LANGUAGES CXX
)

find_package(OpenMP REQUIRED)
find_package(benchmark REQUIRED)
include(CheckCXXSourceRuns)
set(AVX_FLAGS)
set(CMAKE_REQUIRED_FLAGS)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  lib_tests
  src/lib_tests.cpp
)

target_link_libraries(
  lib_tests
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(lib_tests)

set(CMAKE_REQUIRED_FLAGS_SAVED, "${CMAKE_REQUIRED_FLAGS}")

if(MSVC) 
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} /arch:AVX512")
else()
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} -mavx512f -mavx512vl")
endif()

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m512d a, b, c;
          const double src[8] = { 1., 2., 3., 4., 5., 6., 7., 8. };
          double dst[8];
          a =  _mm512_loadu_pd(reinterpret_cast<const __m512d*>(src));
          b =  _mm512_loadu_pd(reinterpret_cast<const __m512d*>(src));
          c = _mm512_add_pd(a, b);
          _mm512_storeu_pd(reinterpret_cast<__m512d*>(dst), c);
          for(int i = 0; i < 8; i++) {
            if((src[i] + src[i]) != dst[i]) {
              return -1;
            }
          }
          return 0;
        }"
        HAVE_AVX512_EXTENSIONS)

set(CMAKE_REQUIRED_FLAGS, "${CMAKE_REQUIRED_FLAGS_SAVED}")

if(MSVC) 
  set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} /O2 /we4309 /we4756 /D_USE_MATH_DEFINES /Wall /wd4711 /wd4710 /wd4127 /wd4577 /wd4100")
else()
  set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -ffast-math -O3 -ftemplate-depth=30000 -ftemplate-backtrace-limit=0 -Werror=overflow -Wall")
endif()

if(HAVE_AVX512_EXTENSIONS)
    message("AVX512 detected, adding benchmarks")
    if(MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512vl")
    endif()
    include_directories(
      benchmarks
      benchmark/include
    )
    add_executable(benchmarks src/benchmarks.cpp src/aerobus.h)
    link_directories(
      benchmarks
      benchmark/build/src
    )
    target_link_libraries(
      benchmarks
      benchmark 
      pthread 
    )
else()
    message("AVX512 not detected, benchmarks need manual configuration")
endif()
