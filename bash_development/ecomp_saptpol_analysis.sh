#!/bin/sh
inputlog=$1
totout="ecomp_summary.dat"

declare -a ecnames=("HarmonicBond" "HarmonicAngle" "CustomTorsion" "PeriodicTorsion" "RBTorsion" "mm.Nonbonded" "Drude" "CustomNonbond")
llength=${#ecnames[@]} #list length

echo $inputlog >> $totout
echo " " >> $totout

i=0
while [ $i -lt $llength ]
do
	ecname=${ecnames[$i]}
	output="temp_${ecname}.dat"

grep $ecname $inputlog | awk '{print $3}' | tail -n +2 > $output
infoline=$(awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {ave = sum[i]/NR; stdev=sqrt((sumsq[i]-sum[i]^2/NR)/NR); 
          printf "%f %f \n", ave, stdev}
         }' $output)
ave="$(cut -d' ' -f1 <<<"$infoline")"
stdev="$(cut -d' ' -f2 <<<"$infoline")"

echo $ecname $ave $stdev >> $totout

	let i=$i+1
done


