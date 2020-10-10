
mv md_npt.pdb md_nvt.pdb
#python ~/des/md_msd.py -i md_nvt.dcd -s start_drudes.pdb -resid_start 0 -resid_end 199 -o chol
python ~/des/md_msd.py -i md_nvt.dcd -s start_drudes.pdb -resid_start 200 -resid_end 399 -o tfn2

