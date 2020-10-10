#!/bin/sh

path=$1
pyscr="/home/kjeong/post_process_codes/hbond_post_matrix_v01.py"
pyexec="/home/kjeong/miniconda3/bin/python"
subdir="test_hbcount"
input="hb_matrix.dat"
hs="histnhb"
sol="nsol"
tstep="2.0"

cd $path
grep 'dcd' *.py
cd $subdir

nohup $pyexec $pyscr -i $input -hs $hs -s $sol -t $tstep &


