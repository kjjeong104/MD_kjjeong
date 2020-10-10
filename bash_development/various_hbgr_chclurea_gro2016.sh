#!/bin/bash

#gr analysis for various pairs
declare -a ht=("cholOHCl" "cholOHuO" "CluHC" "CluHT" "uOuHC" "uOuHT") #hbond type list
#declare -a at1=("uO" "uH" "uN" "cholN" "cholO" "cholHo" "cholHam" "Cl") #atom type list 1
#declare -a at2=() #atom type list 2
declare -a hai1=("15 18" "15 9" "18 19" "18 20" "9 19" "9 20") #atom pair index
#iat1=(9 10 11 13 14 15 17 18) #group index list 1
#iat2=() #group index list 2
llength=${#hai1[@]} #list length

xtcfile=$1
frameeq=$2
frameskip=$3
tprfile="fake_run_forgr_v2016.tpr"
source /home/htjung/.bashrc

#double loop of atom group pair. combination with repetition
#single loop

i=0
while [ $i -lt $llength ]
do
	aname=${ht[$i]}
	outfile="aahbgr_${aname}.xvg"
	cnfile="aahbcn_${aname}.xvg"
	inp="tempindfile"
	echo ${hai1[$i]} | awk '{print $1}'> $inp
	echo ${hai1[$i]} | awk '{print $2}' >> $inp

	gmx rdf -f $xtcfile -b $frameeq -dt $frameskip -cn $cnfile -s $tprfile -n index.ndx -bin 0.01 -excl -o $outfile < $inp

	rm $inp
	i=`expr $i + 1`
done

#while [ $i -lt $llength ]
#do
#	j=$i
#	while [ $j -lt $llength ]
#	do
#		aname1=${ht1[$i]}
#		aname2=${ht1[$j]}
#		outfile="aagr_${aname1}_${aname2}.xvg"
#		cnfile="aacn_${aname1}_${aname2}.xvg"
#		inp="tempindfile"
#		echo ${hai1[$i]} > $inp
#		echo ${hai1[$j]} >> $inp
#		gmx rdf -f $xtcfile -b $frameeq -dt $frameskip -cn $cnfile -s $tprfile -n index.ndx -excl -o $outfile < $inp
#		rm $inp
#		j=`expr $j + 1`
#	done
#
#	i=`expr $i + 1`
#done

