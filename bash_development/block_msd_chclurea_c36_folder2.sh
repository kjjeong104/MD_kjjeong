#!/bin/sh
#xtc file generation

declare -a path=("openmm_c36_400K_m133/" "openmm_c36_400K_m300/" "openmm_c36_400K_m450/")
declare -a ndxs=("../global_reline_m133_c36index.ndx" "../global_reline_m300_c36index.ndx" "../global_reline_m450_c36index.ndx")
declare -a ht=("C2" "cholN" "Cl")
natoms=${#atoms[@]}
declare -a hai1=(12 13 18)
llength=${#hai1[@]}
#path2="nvt_prod"
pyexec="/home/kjeong/miniconda3/bin/python3"
pyscr="/home/kjeong/post_process_codes/dcd_to_xtc_v02.py"

pdbname="nvt_prod.pdb"
dcdname="nvt_prod.dcd"
xtcfile=$1
frameeq=$2
bl=$3 #block length
bn=$4 #block num
#frameskip=$3
#tprfile="../global_fake_run_forgr_c36.tpr"
tprfile="nvt_prod.pdb"
#ndxfile="../global_openmm_c36index.ndx"
source /home/htjung/.bashrc
eftime=`expr $bl / 2`

k=0
for d in "${path[@]}"
do
	cd $d
	#$pyexec $pyscr $pdbname $dcdname $xtcfile
        ndxfile=${ndxs[$k]}
i=0
while [ $i -lt $llength ]
do
	aname=${ht[$i]}
	#outfile="aahbgr_${aname}.xvg"
	#cnfile="aahbcn_${aname}.xvg"
	inp="tempindfile"
	echo ${hai1[$i]} | awk '{print $1}'> $inp
	#echo ${hai1[$i]} | awk '{print $2}' >> $inp

	#gmx msd -f $xtcfile -n $ndxfile -s $tprfile -dt $frameskip -b $frameeq -o msd_${aname}.xvg -endfit $eftime < $inp > log_msd_${aname}.txt
        bindex=0
        while [ $bindex -lt $bn ]
	do
		tb=`expr $bl \* $bindex + $frameeq`
		te=`expr $tb + $bl`
		gmx msd -f $xtcfile -n $ndxfile -s $tprfile -b $tb -e $te -o msd_b${bindex}_${aname}.xvg -endfit $eftime < $inp > log_msd_b${bindex}_${aname}.txt
		let bindex=$bindex+1
	done
        #gmx msd -f $xtcfile -n $ndxfile -s $tprfile -b $frameeq -o msd_${aname}.xvg -endfit $eftime < $inp > log_msd_${aname}.txt
	#gmx rdf -f $xtcfile -b $frameeq -dt $frameskip -cn $cnfile -s $tprfile -n $ndxfile -bin 0.01 -excl -o $outfile < $inp

	rm $inp
	i=`expr $i + 1`
done
	cd ../
        let k=$k+1
done

#for i in "${path1[@]}"
#do
#	#pos="$i/$path2/"
#	pos="$i"
#	pdbname="nvt_prod.pdb"
#	iname="index.ndx"
	#xtc1name="nvt_prod.xtc"
#	xtc2name="nvt_prod_60_110.xtc"
#
#	cd $pos
#	#gmx make_ndx -f $pdbname -o index.ndx < ../../indexgroup_setting.txt
#	j=0
#	while [ $j -lt $natoms ]
#	do
#		atom=${atoms[$j]}
#		iatom=${iat1[$j]}
#		echo $iatom > atomindex_temp.txt
#		gmx msd -f $xtc2name -n $iname -s $pdbname -o msd_${atom}.xvg -endfit $eftime < atomindex_temp.txt > log_msd_${atom}.txt
#		
#                #gmx msd -f $xtc1name -n $iname -b 60000 -e 110000 -s $pdbname -o msd_${atom}.xvg -endfit $eftime < atomindex_temp.txt > log_msd_${atom}.txt
#		rm atomindex_temp.txt
#		let j=$j+1
#	done
	#rm trimset.txt
#	cd ../../
#done
