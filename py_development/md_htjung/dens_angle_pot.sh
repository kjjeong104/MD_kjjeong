# $1: chol/tmpa
# $2: liq/gas
moles=$1
phases=$2
nmol=200

# remove some lines only for potential energy
head -n -3 log > energy.log
sed -i '1,20d' energy.log
sed -i '/class/d' energy.log
awk ' NR % 2 == 0' energy.log > energy_pot.log
awk '{ print $1 }' energy_pot.log > energy_pot.out

# calc mass 
if [ $moles = 'tmpa' ]
then
	mw=$(bc <<< 280.15+102.20)
else 
	mw=$(bc <<< 280.15+104.17)
fi
mass=$(echo $mw*$nmol | bc)

# properties
python ~/des/dens_angle_pot.py $phases $mass | tee dens_angle_pot.log 

