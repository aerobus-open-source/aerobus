set terminal pngcairo size 1200,800
set output 'compensated_horner.png'

set multiplot layout 3,2 title "Evaluation of p(X) = (0.75-X)^5 * (1-X)^{11} (expanded form) in the neighborhood of its multiple roots using Horner and CompHorner\nExact line is the evaluation of the factorized form"  # 2 rows, 3 columns

# Plot 1
set xrange [0.68:1.15]      # Set custom x range
set yrange [-1.5E-12:1.5E-12]     # Set custom y range
set title "Large Sample - Horner"
plot 'large_sample_horner.dat' using 1:2 with lines title 'Horner', \
    'large_sample_horner.dat' using 1:3 with lines title 'Exact'

# Plot 2
set xrange [0.68:1.15]
set yrange [-5E-13:1.5E-12]
set title "Large Sample - CompHorner"
plot 'large_sample_comp_horner.dat' using 1:2 with lines title 'CompHorner', \
    'large_sample_comp_horner.dat' using 1:3 with lines title 'Exact'

# Plot 3
set xrange [0.74995:0.75005]
set yrange [-1.5E-12:1.5E-12]
set title "First Root - Horner"
plot 'first_root_horner.dat' using 1:2 with lines title 'Horner', \
    'first_root_horner.dat' using 1:3 with lines title 'Exact'

# Plot 4
set xrange [0.74995:0.75005]
set yrange [-8.0E-28:8.0E-28]
set title "First Root - CompHorner"
plot 'first_root_comp_horner.dat' using 1:2 with lines title 'CompHorner', \
    'first_root_comp_horner.dat' using 1:3 with lines title 'Exact'

# Plot 5
set xrange [0.9935:1.0065]
set yrange [-1.5E-12:1.5E-12]
set title "Second Root - Horner"
plot 'second_root_horner.dat' using 1:2 with lines title 'Horner', \
    'second_root_horner.dat' using 1:3 with lines title 'Exact'

# Plot 6
set xrange [0.9935:1.0065]
set yrange [-8.0E-28:8.0E-28]
set title "Second Root - CompHorner"
plot 'second_root_comp_horner.dat' using 1:2 with lines title 'CompHorner', \
    'second_root_comp_horner.dat' using 1:3 with lines title 'Exact'

unset multiplot
set output