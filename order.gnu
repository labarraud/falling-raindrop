set terminal postscript eps enhanced color 
set output "error_n2.eps"
set logscale
set xlabel "step h"
set ylabel "error with euclidean norm"
set key on inside right bottom

Eupwind1_n2(x) = ae2*x + be2
fit Eupwind1_n2(x) 'error_upwind1.dat' u (log($1)):(log($2)) via ae2, be2
title_upwind1_n2(ae2, be2) = sprintf('Eupwind1n2(x) = %.2fx + %.2f', ae2, be2)

Eupwind2_n2(x) = ae22*x + be22
fit Eupwind2_n2(x) 'error_upwind2.dat' u (log($1)):(log($2)) via ae22, be22
title_upwind2_n2(ae22, be22) = sprintf('Eupwind2n2(x) = %.2fx + %.2f', ae22, be22)

Eupwind3_n2(x) = ae32*x + be32
fit Eupwind3_n2(x) 'error_upwind3.dat' u (log($1)):(log($2)) via ae32, be32
title_upwind3_n2(ae32, be32) = sprintf('Eupwind3n2(x) = %.2fx + %.2f', ae32, be32)

Ewendroff_n2(x) = aew2*x + bew2
fit Ewendroff_n2(x) 'error_wendroff.dat' u (log($1)):(log($2)) via aew2, bew2
title_wendroff_n2(aew2, bew2) = sprintf('Ewendroffn2(x) = %.2fx + %.2f', aew2, bew2)

Eupwind4_n2(x) = ae42*x + be42
fit Eupwind4_n2(x) 'error_upwind4.dat' u (log($1)):(log($2)) via ae42, be42
title_upwind4_n2(ae42, be42) = sprintf('Eupwind4n2(x) = %.2fx + %.2f', ae42, be42)

plot 'error_upwind1.dat', exp(be2)*x**ae2 t title_upwind1_n2(ae2, be2),'error_upwind2.dat', exp(be22)*x**ae22 t title_upwind2_n2(ae22, be22),'error_upwind3.dat', exp(be32)*x**ae32 t title_upwind3_n2(ae32, be32),'error_upwind4.dat', exp(be42)*x**ae42 t title_upwind4_n2(ae42, be42),'error_wendroff.dat', exp(bew2)*x**aew2 t title_wendroff_n2(aew2, bew2)