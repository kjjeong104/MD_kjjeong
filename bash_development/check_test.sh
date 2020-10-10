#!/bin/sh

i=0
e=999
compstring=" Variable memory released"
joblist="re_joblist.txt"

rm $joblist

while [ $i -le $e ]
do
	out="urea_urea_pbe0_$i.out"
	if [ ! -f "$out" ]; then
		#echo "$i running"
		echo $i >> $joblist
	else
		fin=$(tail -n 1 $out)
		if [ "$fin" == "$compstring" ]; then
			echo "$i complete"
		else
			#echo "$i incomplete terminated"
			echo $i >> $joblist
		fi
	fi
	i=`expr $i + 1`
done
