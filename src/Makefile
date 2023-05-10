clean:
	rm -f *.exe

build:
	g++ -I. -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O0 main.cpp -o main.exe 

run: 
	./main.exe

benchmarks: clean
	g++ -I. -std=c++20 -O3 -mavx2 -mfma -fopenmp benchmarks.cpp -o benchmarks.exe
	./benchmarks.exe

all: clean build run