clean:
	rm -f *.exe

build: # to be sure both compile properly
	clang++ -I. -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 
	g++ -I. -fconcepts-diagnostics-depth=3 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 

debug:
	g++ -I. -std=c++20 -g -ggdb -O0 main.cpp -o main.exe

test: clean build run

run: 
	./main.exe

benchmarks: clean
	g++ -I. -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx512f -mavx512vl -mfma -fopenmp benchmarks.cpp -lm -o benchmarks_avx512.exe
	./benchmarks_avx512.exe

doc:
	rm -rf ../documentation
	cd .. && doxygen Doxyfile
	cd ../documentation/latex && make


all: clean build run benchmarks doc

