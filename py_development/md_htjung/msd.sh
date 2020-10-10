# msd calculation
# $1 : #ion pairs

#mv md_nvt.pdb md_nvt_init.pdb
#mv md_npt.pdb md_nvt.pdb

npairs=$1
id1=$(( $npairs - 1 ))
id2=$(( $id1 + $npairs ))

python ~/des/md_msd.py -i md_nvt_pos.dcd -s start_drudes.pdb -resid_start 0 -resid_end $id1 -o cation
python ~/des/md_msd.py -i md_nvt_pos.dcd -s start_drudes.pdb -resid_start $npairs -resid_end $id2 -o anion

