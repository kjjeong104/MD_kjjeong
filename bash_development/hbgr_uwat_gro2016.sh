#!/bin/bash

#gr analysis for various pairs
declare -a ht=("uO" "uN" "uH" "OW" "HW") #hbond type list
#declare -a at1=("uO" "uH" "uN" "cholN" "cholO" "cholHo" "cholHam" "Cl") #atom type list 1
#declare -a at2=() #atom type list 2
declare -a iat1=(6 7 9 10 11)
#declare -a hai1=("15 18" "15 9" "18 20" "18 21" "9 20" "9 21") #atom pair index
#iat1=(9 10 11 13 14 15 17 18) #group index list 1
#iat2=() #group index list 2
#llength=${#hai1[@]} #list length
llength=${#iat1[@]} #list length
pyexec="/home/kjeong/miniconda3/bin/python3"
pyscr="/home/kjeong/post_process_codes/dcd_to_xtc_v02.py"

pdbname="eq_npt.pdb"
dcdname="eq_npt.dcd"
xtcfile=$1
frameeq=$2
#frameskip=$3
tprfile="fake_run_forgr.tpr"
ndxfile="index.ndx"
source /home/htjung/.bashrc

#loop: get into directory, execute calculations.
#for d in "${path[@]}"
#do
#	cd $d
#	$pyexec $pyscr $pdbname $dcdname $xtcfile
i=0
while [ $i -lt $llength ]
do
	j=$i
	while [ $j -lt $llength ]
	do
		aname1=${ht[$i]}
		aname2=${ht[$j]}
		outfile="aahbgr_${aname1}_${aname2}.xvg"
		cnfile="aahbcn_${aname1}_${aname2}.xvg"
		inp="tempindfile"
		echo ${iat1[$i]} > $inp
		echo ${iat1[$j]} >> $inp

		gmx rdf -f $xtcfile -b $frameeq -cn $cnfile -s $tprfile -n $ndxfile -bin 0.01 -excl -o $outfile < $inp

		rm $inp
		j=`expr $j + 1`
	done
	i=`expr $i + 1`
done
#	cd ../
#done

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

