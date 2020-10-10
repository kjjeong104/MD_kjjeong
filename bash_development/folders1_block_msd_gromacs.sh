#!/bin/sh

declare -a path1=("descale030_ge1p_moexp80hel80/" "size72_nvt_m030_m80h80/" "size288_nvt_m030_m80h80/")
declare -a atoms=("CU" "OW")
natoms=${#atoms[@]}
iat1=(8 10)
source /home/htjung/.bashrc
teq=5000 #eqtime : 5ns
bt=2000 #block time: 2ns
nb=5 #number of blocks:5
tfit=5000 #total MSD fit : cut until 5ns
sfit=2000 #shorter fit
tbfit=1000 #block MSD fit : cut until 1ns

for path in "${path1[@]}"
do
	pdbname="eq_nvt.pdb"
	iname="index.ndx"
	xtcname="eq_nvt.xtc"

	cd $path
	tprname=$(ls *.tpr)
	j=0
	while [ $j -lt $natoms ]
	do
		atom=${atoms[$j]}
		iatom=${iat1[$j]}
		#total MSD: cut half.
		echo $iatom | gmx msd -f $xtcname -n $iname -s $tprname -b $teq -endfit $tfit -o msd_tfit_${atom}.xvg  > log_msdtfit_${atom}.txt
		echo $iatom | gmx msd -f $xtcname -n $iname -s $tprname -b $teq -endfit $sfit -o msd_sfit_${atom}.xvg  > log_msdsfit_${atom}.txt

		k=0
		while [ $k -lt $nb ]
		do
			bbt=`expr $teq + $k \* $bt`
			bet=`expr $teq + $k \* $bt + $bt`
			echo $iatom | gmx msd -f $xtcname -n $iname -s $tprname -b $bbt -e $bet -endfit $tbfit -o msd_${atom}_b${k}.xvg >> log_bmsdfit_${atom}.txt

			let k=$k+1
		done

		let j=$j+1
	done
	cd ../
done
