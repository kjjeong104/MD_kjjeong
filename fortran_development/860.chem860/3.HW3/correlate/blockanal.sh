#!/bin/sh

nblock=20
destin="../mysol1anew"
i=1
while [ $i -le $nblock ]
do
	cp ../vel_$i.out ./velocity.dat
	./correlate_original < quick.inp
	
	mv acf.out $destin/acf_$i.out
	mv corr.out $destin/corr_$i.out
	i=`expr $i + 1`
done
