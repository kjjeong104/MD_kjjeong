#!/bin/sh

nptlog=$1
vollog="npt_volume.log"
frameskip=$2
let frameskip=$frameskip+1
avefile="ave_volume.log"

#grep ' nm' $nptlog | awk '{if(NR%2==0) print $1}' > $vollog
grep ' nm' $nptlog | grep -v 'kJ' > $vollog
avevol=$(tail -n +$frameskip $vollog | awk '{sum+=$1}END{print sum/NR}')
boxsize=$(awk -v num=$avevol 'BEGIN{print num^(1/3)}')

echo $avevol >> $avefile
echo " " >> $avefile
echo $boxsize >> $avefile
