#!/bin/sh
#script to write input file for cpptraj, 

input="combine_crds.in"
topfile="oplsdes_chclurea_m200.prmtop"
prefix=$1 #mdcrd file prefix
num=$2 #number of files

if [ -f $input ];then
	rm $input
fi

i=1
while [ $i -le $num ]
do
	if [ $i -eq 1 ];then
		echo "trajin ${prefix}.1GPU" >> $input
	else
		echo "trajin ${prefix}${i}.1GPU" >> $input
	fi
	let i=$i+1
done
echo "trajout combine_${prefix}.xtc" >> $input
echo "go" >> $input

cpptraj -p $topfile < $input

