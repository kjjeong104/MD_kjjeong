mkdir 2nvt
cd 2nvt
cp ../../323k/2nvt/run_openMM_nvt.py ./
cp ../1nvt/md_nvt.pdb ./
sed -i 's/TMP /TMPA/g' md_nvt.pdb
sed -i 's/Tf2 /Tf2N/g' md_nvt.pdb
sed -i '/ D/d' md_nvt.pdb
sed -i '/CONECT/d' md_nvt.pdb
cp ../1nvt/md_nvt.chk ./
cp ~/des/sapt* ./

