### this makefile is used for development -- use cmake

clean:
	rm -f *.exe
	rm -rf documentation
	rm -rf build

doc:
	doxygen Doxyfile
	cd documentation/latex && make

tests:
	cmake -S . -B build
	cmake --build build
	cd build && ctest

all: clean build run

