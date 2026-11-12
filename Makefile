### this makefile is used for development -- use cmake

.PHONY: clean doc lint check_test_names tests all paper

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

# a test file is tests/<topic>_tests.cpp; the other directories under tests/ hold their own kinds of
# sources (compile_fail/, cuda/, hip/, hygiene/, package/) and are left alone
check_test_names:
	@bad=$$(find tests -maxdepth 1 -name '*.cpp' ! -name '*_tests.cpp'); \
	if [ -n "$$bad" ]; then \
		echo "test files must be named tests/<topic>_tests.cpp:"; \
		echo "$$bad" | sed 's/^/  /'; \
		exit 1; \
	fi

# the same files the cpplint job of the CI lints
lint: check_test_names
	cpplint src/aerobus.h src/benchmarks.cpp tests/*.cpp examples/*.cpp

tests: clean
	echo "RUNNING TESTS WITH CLANG"
	cmake -D CMAKE_C_COMPILER=clang -D CMAKE_CXX_COMPILER=clang++ -S . -B build
	cmake --build build
	cd build && ctest
	rm -rf build
	echo "RUNNING TESTS WITH GCC"
	cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -S . -B build 
	cmake --build build
	cd build && ctest

all: lint tests doc

paper: clean
	cd paper && pdflatex main.tex
	cd paper && bibtex main
	cd paper && pdflatex main.tex

