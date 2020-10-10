#!/bin/sh

fnum=1
while [ -f "test${fnum}.txt" ]
do
	let fnum=$fnum+1
done
echo $fnum
