#!/bin/sh
mdlog=$1
temp='temp_energy.txt'
frameskip=$2
let frameskip=$frameskip+1
#corrtxt='tcorr_PE.txt'
#inve=0.367879

#step 1: average and stdev of potential energy
grep 'kJ/mol' $mdlog | grep -v '<' | awk 'NR%2==0' | awk '{print $1}' | tail -n +3 | tail -n +$frameskip > $temp
#grep 'kJ/mol' $mdlog | awk 'NR % 10 == 2' | awk '{print $1}' | tail -n +3 | tail -n +$frameskip > $temp
echo "pot energy : ave and stdev, in kJ/mol"
infoline=$(awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}} 
          END {for (i=1;i<=NF;i++) {ave = sum[i]/NR; stdev=sqrt((sumsq[i]-sum[i]^2/NR)/NR);
          printf "%f %f \n", ave, stdev}
         }' temp_energy.txt)
ave="$(cut -d' ' -f1 <<<"$infoline")"
stdev="$(cut -d' ' -f2 <<<"$infoline")"
echo $ave $stdev

#rm $temp

