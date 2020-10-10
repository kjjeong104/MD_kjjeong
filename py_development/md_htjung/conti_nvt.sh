mkdir 2nvt
cd 1nvt/
cp sapt* ../2nvt/
cp *.chk ../2nvt/
cp nvt_init.pdb ../2nvt/
cp run_openMM_nvt.py ../2nvt/
cd ..
cd 2nvt
sed -i "s/#simmd.loadCheckpoint/simmd.loadCheckpoint/g" run_openMM_nvt.py
sed -i "s/simmd.step(100000)/#simmd.step(100000)/g" run_openMM_nvt.py
sed -i "/velocity stable/d" run_openMM_nvt.py
sed -i "s/range(1,10001)/range(1,20001)/g" run_openMM_nvt.py

