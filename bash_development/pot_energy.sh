#!/bin/sh
mdlog='md.log'
temp='temp_energy.txt'
corrtxt='tcorr_PE.txt'
inve=0.367879

#step 1: average and stdev of potential energy
grep 'kJ/mol' $mdlog | awk 'NR % 9 == 2' | awk '{print $1}' | tail -n +3 > $temp
echo "pot energy : ave and stdev, in kJ/mol"
infoline=$(awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {ave = sum[i]/NR; stdev=sqrt((sumsq[i]-sum[i]^2/NR)/NR);
          printf "%f %f \n", ave, stdev}
         }' temp_energy.txt)
ave="$(cut -d' ' -f1 <<<"$infoline")"
stdev="$(cut -d' ' -f2 <<<"$infoline")"
echo $ave $stdev

#step 2: calculating autocorrelation
nline=$(wc -l $temp | awk '{print $1}')
trange=`expr $nline \/ 2`
declare -a array1
declare -a ct #autocorr fxn array
declare -a ctcount
while read line
do
	array1+=("$line")
done < $temp
i=0
ct[0]=1.000
ctcount[0]=1.000
while [ $i -lt $trange ]
do
	k=1
	echo $i
	while [ $k -lt $trange ]
	do
		j=`expr $i + $k`
		ct[k]=$(awk -v c="${ct[$k]}" -v x1="${array1[$i]}" -v x2="${array1[$j]}" -v m="$ave" -v s="$stdev" 'BEGIN {print c + (x1-m)(x2-m)/(s*s)}')
		ctcount[k]=$(awk -v count="${ctcount[$k]}" 'BEGIN {print count + 1}' )
		k=`expr $k + 1`
	done
	i=`expr $i + 1`
done

i=0
e=${#ct[@]}
while [ $i -lt $e ]
do
	ct[i]=$(awk -v a="${ct[$i]}" -v b="${ctcount[$i]}" 'BEGNI {print a / b}' )
	i=`expr $i + 1`
done

#echo "rough corrt in steps"
#echo $tau

#print time corr fxn into output
printf '%s\n' "${ct[@]}" 

rm $temp
