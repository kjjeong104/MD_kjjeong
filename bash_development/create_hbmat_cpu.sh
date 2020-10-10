#!/bin/sh

path=$1
setfile="../../settings_hbmat_NHO.txt"
pyscr="/home/kjeong/post_process_codes/hbond_matrix_v01.py"
pyexec="/home/kjeong/miniconda3/bin/python"
subdir="test_hbcount"
outfile="hb_matrix.dat"
logfile="log_hbmat_1.txt"
pdbfile="../eq_npt.pdb"
dcdfile="../eq_npt.dcd"

cd $path
if [ -d $subdir ]
then
	echo "accessible"
else
	mkdir $subdir
fi
cd $subdir

nohup $pyexec $pyscr $pdbfile $dcdfile $outfile < $setfile > $logfile &


