set terminal pngcairo size 1200,800
set output 'double_vs_comp_float.png'

set multiplot layout 1,2 title "Evaluation of p(X) = (x-1)^{3} (expanded form) in the neighborhood of its root using Horner (double) and CompHorner (float)\nExact line is the evaluation of the factorized form"  # 2 rows, 3 columns

# Plot 1
set xrange [0.99995:1.00005]
set yrange [-2E-13:2E-13]
set title "Horner Double"
plot 'double.dat' using 1:2 with lines title 'Double', \
    'double.dat' using 1:3 with lines title 'Exact'

# Plot 2
set xrange [0.99995:1.00005]
set yrange [-2E-13:2E-13]
set title "Comp Float"
plot 'float_comp.dat' using 1:2 with lines title 'CompFloat', \
    'float_comp.dat' using 1:3 with lines title 'Exact'

unset multiplot
set output