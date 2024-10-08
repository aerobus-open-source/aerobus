rm -f plots/*.dat plots/*.png a.out
clang++ -std=c++20 -O3 -mfma -mavx512f -mavx512vl compensated_horner.cpp
./a.out
cd plots && gnuplot compensated_horner.gp && gnuplot double_horner_vs_float_comp_horner.gp