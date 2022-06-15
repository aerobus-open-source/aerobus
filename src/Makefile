clean:
	rm -f main.exe

build:
	clang++ -I. -std=c++20 -O3 -mavx2 -mfma -fopenmp main.cpp -o main.exe

run: 
	./main.exe

all: clean build run