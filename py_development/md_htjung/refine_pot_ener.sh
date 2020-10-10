# refine pot energy from log file

cp md*.log energy_pot.log

sed -i '/CUDA/d' energy_pot.log
sed -i '/Method/d' energy_pot.log
sed -i '/openmm/d' energy_pot.log
sed -i '/TFSI/d' energy_pot.log
sed -i '/NVT/d' energy_pot.log
sed -i '/Simu/d' energy_pot.log
sed -i '/sim/d' energy_pot.log
sed -i '/Done/d' energy_pot.log

awk 'NR % 2 == 0' energy_pot.log > energy_pot.bak
mv energy_pot.bak energy_pot.log

sed -i -e 's@kJ/mol@ @g' energy_pot.log

