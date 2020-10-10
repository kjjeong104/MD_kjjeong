#!/bin/sh
#Shell script for ensemble-averaging the dE/dlambda

input=$1
nstep=500
nlblock=11 #number of lambda blocks

typestr=$(head -n 1 $input | awk '{print $7}')
output="ensdedl_${typestr}.dat"

i=0
while [ $i -lt $nlblock ]
do
	#lambda line n1: 2+($nstep+2)*i. Block start n2: n1+2, Block end n3: n2+nstep-1
	let n1="$i*($nstep+2)+2"
	let n2=$n1+2
	let n3=$n2+$nstep-1
	temp="blockenergy_temp.dat"
	sed -n "${n2},${n3}p" $input | awk '{print $3}' > $temp
	lambda=$(sed -n "${n1},${n1}p" $input | awk '{print $5}')
	infoline=$(awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {ave = sum[i]/NR; stdev=sqrt((sumsq[i]-sum[i]^2/NR)/NR); 
          printf "%f %f \n", ave, stdev}
         }' $temp)
	dedlave="$(cut -d' ' -f1 <<<"$infoline")"
	dedlstdev="$(cut -d' ' -f2 <<<"$infoline")"

	echo "$lambda $dedlave $dedlstdev" >> $output
	rm $temp

	let i=$i+1
done

