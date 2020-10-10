#!/bin/sh

input=$1

i=1
end=2
period=306

while [ $i -le $end ]
do
	output="snapshot_FPMD_333_AnR_${i}.pdb"
	j1=`expr $i \* $period - $period + 1`
	j2=`expr $i \* $period`

	sed -n "${j1},$j2" $input > $output
	let i=$i+1
    
        if (( $i % 100 == 0))then
		echo "snapshot $i"
	fi
done

