#!/bin/bash
# cat dcd files
# $1: initial number for folder
# $2: final number for folder
# $3: prefix of output dcd file name

prefix_folder="eq"
output=$3
init=$1
final=$2
file1="md_nvt"
log1="md.log"
dcd=".dcd"
space=" "
comline=""
slash="/"

# cat trajectory
start=$init
until [ $start -gt $final ]
do
	temp=$prefix_folder$start$slash$file1$dcd
	comline=$comline$temp$space
	let start=start+1
done
~/catdcd4.0/catdcd -o $output$dcd $comline

# cat log files
start=$init
until [ $start -gt $final ]
do
	cat $prefix_folder$start$slash$log1 >> md-cat.log
	let start=start+1
done

