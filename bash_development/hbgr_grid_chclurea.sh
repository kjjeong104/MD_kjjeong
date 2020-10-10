#!/bin/bash

#gr analysis for various pairs
#declare -a path=("./")
declare -a ht=("cholOHCl" "cholOHuO" "CluHC" "CluHT" "uOuHC" "uOuHT") #hbond type list
#declare -a at1=("uO" "uH" "uN" "cholN" "cholO" "cholHo" "cholHam" "Cl") #atom type list 1
#declare -a at2=() #atom type list 2
declare -a hai1=("15 18" "15 9" "18 20" "18 21" "9 20" "9 21") #atom pair index
#iat1=(9 10 11 13 14 15 17 18) #group index list 1
#iat2=() #group index list 2
llength=${#hai1[@]} #list length
pyexec="/home/kjeong/miniconda3/bin/python3"
pyscr="/home/kjeong/post_process_codes/dcd_to_xtc_v02.py"
gridsize=1000
gname="grid"
pdbname="prod303_nvt.pdb"
dcdname="prod303_nvt.dcd"
xtcfile="prod303_nvt.xtc"
frameeq=1000
#frameskip=$3
tprfile="../global_fake_ureaom_forgr.tpr"
ndxfile="../global_ureaom_index.ndx"
source /home/htjung/.bashrc

g=0
#loop: get into directory, execute calculations.
while [ $g -lt $gridsize ]
#for d in "${path[@]}"
do
	d="$gname$g"
	cd $d
	$pyexec $pyscr $pdbname $dcdname $xtcfile
i=0
while [ $i -lt $llength ]
do
	aname=${ht[$i]}
	outfile="aahbgr_${aname}.xvg"
	cnfile="aahbcn_${aname}.xvg"
	inp="tempindfile"
	echo ${hai1[$i]} | awk '{print $1}'> $inp
	echo ${hai1[$i]} | awk '{print $2}' >> $inp

	gmx rdf -f $xtcfile -b $frameeq -cn $cnfile -s $tprfile -n $ndxfile -bin 0.01 -excl -o $outfile < $inp

	rm $inp
	i=`expr $i + 1`
done
	cd ../
	g=`expr $g + 1`
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

