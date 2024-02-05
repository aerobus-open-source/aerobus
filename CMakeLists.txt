cmake_minimum_required(VERSION 3.27)

project(
    Aerobus
    VERSION 1.2
    LANGUAGES CXX
)

find_package(OpenMP REQUIRED)

set(AVX_FLAGS)
set(CMAKE_REQUIRED_FLAGS)
include(CheckCXXSourceRuns)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m256 a = _mm256_fmadd_ps(
            _mm256_set1_ps(0.0F),
            _mm256_set1_ps(0.0F),
            _mm256_set1_ps(0.0F));
          return 0;
        }"
        HAVE_AVX_EXTENSIONS)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512vl")

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m512 a = _mm512_fmadd_ps(
            _mm512_set1_ps(0.0F),
            _mm512_set1_ps(0.0F),
            _mm512_set1_ps(0.0F));
          return 0;
        }"
        HAVE_AVX512_EXTENSIONS)


set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if(HAVE_AVX512_EXTENSIONS)
    message("AVX512 detected, adding flags")
    if(MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512vl -mfma")
    endif()
elseif(HAVE_AVXFMA_EXTENSIONS)
    message("AVX2+FMA detected, adding flags")
    if(MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2 /fp:fast")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma")
    endif()
else() 
    message("not vectorization detected - skip benchmarks")
endif()


if(MSVC) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D__GALOIS_FIELDS__ /Ot /Ob2 /O2 /we4309 /we4756 /D_USE_MATH_DEFINES /Wall /wd4711 /wd4710 /wd4127 /wd4577 /wd4100")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -D__GALOIS_FIELDS__ -Wno-unused-local-typedefs -ftemplate-depth=30000 -ftemplate-backtrace-limit=0 -Werror=overflow -Wall")
endif()

add_executable(aerobus_test src/main.cpp src/lib.h imports/conwaypolynomials.h)

if(HAVE_AVXFMA_EXTENSIONS OR HAVE_AVX512_EXTENSIONS)
    add_executable(aerobus_benchmarks src/benchmarks.cpp src/lib.h imports/conwaypolynomials.h)
endif()