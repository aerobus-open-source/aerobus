cmake_minimum_required(VERSION 3.1...3.27)

project(
    Aerobus,
    VERSION 1.2
    LANGUAGES CXX
)

find_package(OpenMP REQUIRED)
include(CheckCXXSourceRuns)
set(AVX_FLAGS)
set(CMAKE_REQUIRED_FLAGS)

check_cxx_source_runs("
        #include <immintrin.h>
        int main()
        {
          __m512d a, b, c;
          const double src[8] = { 1., 2., 3., 4., 5., 6., 7., 8. };
          double dst[8];
          a =  _mm512_loadu_pd( (__m512d*)src );
          b =  _mm512_loadu_pd( (__m512d*)src );
          c = _mm512_add_pd( a, b );
          _mm512_storeu_pd( (__m512d*)dst, c );
          for( int i = 0; i < 8; i++ ){
            if( ( src[i] + src[i] ) != dst[i] ){
              return -1;
            }
          }
          return 0;
        }"
        HAVE_AVX512_EXTENSIONS)

set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -O3 -ftemplate-depth=30000 -ftemplate-backtrace-limit=0 -Werror=overflow -Wall")
set(CMAKE_CXX_STANDARD 20)

add_executable(aerobus_test src/main.cpp)
if(HAVE_AVX512_EXTENSIONS)
    message("AVX512 detected, adding benchmarks")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512f -mavx512vl")
    add_executable(aerobus_benchmarks src/benchmarks.cpp)
else()
    message("AVX512 not detected, benchmarks need manual configuration")
endif()
