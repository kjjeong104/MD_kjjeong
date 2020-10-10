# extract pdb files from dcd files every steps

read_name="md_nvt.dcd"
top_name="start_drudes.pdb"
time_step=100
init_time=5000
n_pdbs=10

counter=0
iframe=0
while [ $counter -lt $n_pdbs ];do
	iframe=$(( $counter * $time_step + $init_time))
	echo $iframe
	python ~/des/trj_to_pdb.py $read_name $iframe $top_name
	sed -i 's/Cho /Chol/g' out.pdb
	sed -i 's/TMP /TMPA/g' out.pdb
	sed -i 's/Tf2 /Tf2N/g' out.pdb
	sed -i '/         D/d' out.pdb
	sed -i '/CONECT/d' out.pdb
        mv out.pdb nvt_$counter.pdb
	let counter=counter+1
done

