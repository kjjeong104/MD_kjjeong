python ~/des/dcd_pos_vel_combine.py md_nvt_pos.dcd md_nvt_vel.dcd start_drudes.pdb
gmx editconf -f conf.gro -resnr 1 -o conf.gro
~/des/msd.sh
python ~/des/viscosity_helfand.py -i traj.trr -s conf.gro -temp $1 -step 5

