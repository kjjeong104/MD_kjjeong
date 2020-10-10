#!/bin/bash
# convert dcd file to trr and cat trajectory with gromacs
# $1: initial number for folder
# $2: final number for folder
# you should install mdtraj

prefix_folder="eq"
output="md"
init=$1
final=$2
file1="md_npt"
dcd=".dcd"
trr=".trr"
space=" "
slush="/"
down="../"
gmx="gmx_mpi"
hypen="-"

# convert dcd to trr
start=$init
until [ $start -gt $final ]
do
	move=$prefix_folder$start
	cd $move
	mdconvert -o $down$move$trr $file1$dcd
	cd ..
	let start=start+1
done
 
# cat .trr
$gmx trjcat -f *.trr -settime -o $output$init$hypen$final.trr
rm $prefix_folder*.trr

