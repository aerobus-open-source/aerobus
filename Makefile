### this makefile is used for development -- use cmake

clean:
	rm -f *.exe
	rm -rf documentation
	rm -rf build
	rm -f paper/*.aux paper/*.bbl paper/*.blg paper/*.log paper/*.pdf paper/*.out

# version from src/aerobus.h, used by the Doxyfile as $(AEROBUS_VERSION)
AEROBUS_VERSION := v$(shell sed -n 's/^\#define AEROBUS_VERSION_\(MAJOR\|MINOR\|PATCH\) //p' src/aerobus.h | paste -sd.)

doc:
	AEROBUS_VERSION=$(AEROBUS_VERSION) doxygen Doxyfile
	rm -rf docs/ && mkdir -p docs && cp -r documentation/html/* docs/

tests:
	echo "RUNNING TESTS WITH CLANG"
	cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -S . -B build
	cmake --build build
	cd build && ctest
	rm -rf build
	echo "RUNNING TESTS WITH GCC"
	cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -S . -B build 
	cmake --build build
	cd build && ctest

all: tests doc

paper: clean
	cd paper && pdflatex main.tex
	cd paper && bibtex main
	cd paper && pdflatex main.tex

