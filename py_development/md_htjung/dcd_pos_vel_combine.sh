# dcd to trr
#  then calc viscosity

nstart=$1
nend=$2

while [ $nstart -lt $nend ];do
	mv md_nvt_pos.dcd.$nstart md_nvt_pos.dcd
	mv md_nvt_vel.dcd.$nstart md_nvt_vel.dcd
	mv start_drudes.pdb.$nstart start_drudes.pdb
	python ~/des/dcd_pos_vel_combine.py md_nvt_pos.dcd md_nvt_vel.dcd start_drudes.pdb
	echo 0 | gmx tcaf -f traj.trr -s conf.gro -ov >> vis_log
	rm conf.gro traj.trr
	let nstart=nstart+1
done

