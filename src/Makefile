clean:
	rm -f *.exe

build: # to be sure both compile prolerly
	clang++ -I. -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 
	g++ -I. -fconcepts-diagnostics-depth=3 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 

run: 
	./main.exe

benchmarks: clean
	g++ -I.  -DDEGREE=1 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=3 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=5 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=7 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=9 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=11 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=13 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=15 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe
	g++ -I.  -DDEGREE=17 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe

all: clean build run