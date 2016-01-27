set terminal postscript eps enhanced color 
set output "plotVxvelocity.eps"
set pm3d map
splot 'particleinit.dat' matrix