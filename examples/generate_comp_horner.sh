clang++ -std=c++20 compensated_horner.cpp
rm -f plots/*.dat plots/*.png
./a.out
cd plots && gnuplot compensated_horner.gp