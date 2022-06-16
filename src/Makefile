clean:
	rm -f *.exe

build:
	clang++ -I. -std=c++20 -O3 -mavx2 -mfma -fopenmp main.cpp -o main.exe

run: 
	./main.exe

benchmarks: clean
	clang++ -I. -std=c++20 -O3 -mavx2 -mfma -fopenmp benchmarks.cpp -o benchmarks.exe
	./benchmarks.exe

all: clean build run