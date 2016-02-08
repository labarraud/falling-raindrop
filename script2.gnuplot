set terminal postscript eps enhanced color 
set output "plotparticleinit.eps"
set pm3d map
splot 'particleinit.dat' matrix

set terminal postscript eps enhanced color 
set output "plotparticlefinal.eps"
set pm3d map
splot 'particlefinal.dat' matrix
