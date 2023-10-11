clean:
	rm -f *.exe

build:
	g++ -I. -fconcepts-diagnostics-depth=3 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 

run: 
	./main.exe

benchmarks: clean
	g++ -I../../VML -I. -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 -mavx2 -mfma -fopenmp benchmarks.cpp -lm -o benchmarks.exe
	./benchmarks.exe

all: clean build run