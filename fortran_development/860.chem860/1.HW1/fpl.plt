set size 1,1
set nogrid
set nokey
set nopolar
set angles radians

set terminal x11
set autoscale y
set autoscale x
set title "" 0,-8
set xlabel 'r' 0,.5
set xtics ('0' 0, '-.2' -.2, '.2' .2, '-.4' -.4, '.4' .4)
set ylabel 'p'  0.5,0
set zero 1e-08

#################################
plot 'nm_time.dat' using 2:3 wi li lt 1 lw 3
pause -1
set terminal postscript
set output 'fpl.ps'
replot
set output 'STDOUT'



quit
