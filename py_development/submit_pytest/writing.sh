#!/bin/sh

# untar your Python installation
tar -xzf python.tar.gz
# make sure the script will use your Python installation, 
# and the working directory as it's home location
export PATH=$(pwd)/python/bin:$PATH
mkdir home
export HOME=$(pwd)/home
# run your script
#python my_script.py

num=$1

initt=100
bsize=50

resultt=`expr $initt + $bsize \* $num`

infile=input_$resultt.txt
outfile=output_$resultt.txt

python memo10.py $infile $outfile
