#!/bin/bash
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --time=3-00:00
##SBATCH --mem=6000M
##SBATCH --partition=gpudefault
#SBATCH --partition=yethiraj
##SBATCH --constraint=1080ti
##SBATCH --constraint=gtx2070
#SBATCH --job-name=test
#SBATCH --output=slurm_test1
#SBATCH --error=slurm_test2

pyexec="/home/kjeong/miniconda3/bin/python"
pyscr="template_run_apt1t3.py"
log="md_apt1t3.log"
dirindex=$1

source ~/.bashrc
startdir="/home/kjeong/3.Deep_Eutectic_Solvent/c1.ChCl_urea_12/expden_303K_runs/gridtest_saptff"
dir="$startdir/grid${dirindex}"

cd $dir

$pyexec $pyscr > $log

