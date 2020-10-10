#!/bin/sh

input="dihscan_total.xyz"
pdbform="choline.pdb"
temppdb="temp.pdb"
lineran=20
set1=24
ang=0
i=0
inter=3
maxang=360
temp="temp.xyz"

sed -n "2,22p" $pdbform > $temppdb

while [ $ang -lt $maxang ]
do
	output="qm_HCCO_${ang}.pdb"
	l1=`expr $i \* $set1 + 3`
	l2=`expr $l1 + $lineran`
	head -n 1 $pdbform >> $output
	sed -n "$l1,${l2}p" $input | awk '{print $2 " " $3 " " $4}' > $temp
	
	awk '
		NR==FNR {
		A[FNR]=sprintf("%7.3f %7.3f %7.3f", $1, $2, $3)
		next
		}
		{
		print substr($0,1,p) A[FNR] substr($0,p+length(A[FNR])+1)
		}
	' p=31 $temp $temppdb >> $output
	tail -n 1 $pdbform >> $output

	ang=`expr $ang + $inter`
	i=`expr $i + 1`
done

rm $temp $temppdb
