set terminal postscript eps enhanced color 
set output "plotvelocity.eps"
set xlabel "x"
set ylabel "y"
set xrange [0:5]
set yrange [0:5]
set key off
plot 'velocity.dat' w vec