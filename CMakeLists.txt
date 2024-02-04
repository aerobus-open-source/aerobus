cmake_minimum_required(VERSION 3.1...3.27)

project(
    Aerobus
    VERSION 1.2
    LANGUAGES CXX
)

find_package(OpenMP REQUIRED)

include(CheckCXXSourceRuns)

set(CMAKE_REQUIRED_FLAGS)

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m512d a, b, c;
          const double src[8] = { 1., 2., 3., 4., 5., 6., 7., 8. };
          double dst[8];
          a =  _mm512_loadu_pd(src);
          b =  _mm512_loadu_pd(src);
          c =  _mm512_loadu_pd(src);
          d = _mm512_fmadd_pd(a, b, c);
          _mm512_storeu_pd(dst, c);
          for(int i = 0; i < 8; i++) {
            if((src[i] * src[i] + src[i]) != dst[i]) {
              return -1;
            }
          }
          return 0;
        }"
        HAVE_AVX512_EXTENSIONS)

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m256d a, b, c, d;
          const double src[4] = { 1., 2., 3., 4. };
          double dst[4];
          a =  _mm256_loadu_pd(src);
          b =  _mm256_loadu_pd(src);
          c =  _mm256_loadu_pd(src);
          d = _mm256_fmadd_pd(a, b, c);
          _mm256_storeu_pd(dst, d);
          for(int i = 0; i < 4; i++) {
            if((src[i] * src[i] + src[i]) != dst[i]) {
              return -1;
            }
          }
          return 0;
        }"
        HAVE_AVX2_FMA_EXTENSIONS)


set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

if(HAVE_AVX512_EXTENSIONS)
    message("AVX512 detected, adding flags")
    if(MSVC)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512vl -mfma")
    endif()
elseif(HAVE_AVX2_FMA_EXTENSIONS)
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
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -D__GALOIS_FIELDS__ -ftemplate-depth=30000 -ftemplate-backtrace-limit=0 -Werror=overflow -Wall")
endif()

message("CXX_FLAGS : ${CMAKE_CXX_FLAGS}")

add_executable(aerobus_test src/main.cpp src/lib.h imports/conwaypolynomials.h)

if(HAVE_AVX2_FMA_EXTENSIONS OR HAVE_AVX512_EXTENSIONS)
    add_executable(aerobus_benchmarks src/benchmarks.cpp src/lib.h imports/conwaypolynomials.h)
endif()