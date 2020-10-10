# npt to nvt simulation prepration
# $1: chol/tmpa/line number from density.log

read_name="md_npt.dcd"
top_name="md_npt.pdb"
moln=$1
ofolder="../2nvt/"

# get density info or pick iframe
if [ "$moln" == "chol" ]; then 
	~/des/dens_angle_pot.sh $moln liq
	pdb_idx=$(grep 'index' dens_angle_pot.log | awk '{ print $7 }')
	let pdb_idx=pdb_idx+1
elif [ "$moln" == "tmpa" ]; then
	~/des/dens_angle_pot.sh $moln liq
	pdb_idx=$(grep 'index' dens_angle_pot.log | awk '{ print $7 }')
	let pdb_idx=pdb_idx+1
else
	pdb_idx=$moln	
fi
python ~/des/trj_to_pdb.py $read_name $pdb_idx $top_name
# out.pdb file is what we want to run in NVT
mv out.pdb nvt_init.pdb
sed -i 's/Cho /Chol/g' nvt_init.pdb
sed -i 's/Tf2 /Tf2N/g' nvt_init.pdb
sed -i 's/TMP /TMPA/g' nvt_init.pdb
sed -i '/         D/d' nvt_init.pdb
sed -i '/CONECT/d' nvt_init.pdb
mkdir $ofolder
cp nvt_init.pdb $ofolder
cp *.xml $ofolder
cp ~/des/run_openMM_nvt.py $ofolder

