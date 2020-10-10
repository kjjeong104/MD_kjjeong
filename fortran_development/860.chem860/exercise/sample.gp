set size 1,1	#size of plot
set nogrid
set nokey
set nopolar
set angles radians

#
#
set terminal x11
set output 'STDOUT'

#
#

#
#

set autoscale y
set autoscale x

#
#
#

set title "" 0,-8 #title
set xlabel 'x' 0, .5
set ylabel 'y' 0.5,0
set zero 1e-08

###
#
plot 'file_name.dat' using 1:2 wi imp lt 1 lw 1, 'file_name2.dat' using 1:2 wi li lt 2 lw 2
pause -1

quit

