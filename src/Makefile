clean:
	rm -f *.exe

build-profile:
	../../llvm-project/build/bin/templight++ -Xtemplight --profiler -Xtemplight -output=main -Xtemplight -blacklist=templight.blacklist -Xtemplight -ignore-system -I. -ftemplate-depth=10000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 
	../../templight-tools/build/bin/templight-convert --format callgrind --output trace.cg main.trace.pbf

build: # to be sure both compile properly
	clang++ -D__GALOIS_FIELDS__ -I. -ftemplate-depth=10000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 
	# g++ -D__GALOIS_FIELDS__ -I. -fconcepts-diagnostics-depth=3 -ftemplate-depth=30000 -std=c++20 -ftemplate-backtrace-limit=0 -O3 main.cpp -o main.exe 

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

