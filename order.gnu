set terminal postscript eps enhanced color 
set output "error_n2.eps"
set logscale
set xlabel "step h"
set ylabel "error with euclidean norm"
set key on inside right bottom

Eeuler_n2(x) = ae2*x + be2
fit Eeuler_n2(x) 'error_upwind1.dat' u (log($1)):(log($2)) via ae2, be2
title_Eeuler_n2(ae2, be2) = sprintf('Eeulern2(x) = %.2fx + %.2f', ae2, be2)

plot 'error_upwind1.dat', exp(be2)*x**ae2 t title_Eeuler_n2(ae2, be2)