# $1 : coh/no

if [ $1 = 'coh' ]; then
	grep "CustomNonbondedForce" log | awk '{ print $3 }' > pot_cnbond.log
	grep "NonbondedForce" log | awk ' NR % 2 == 1' | awk '{ print $3 }' > pot_nbond.log
else
	grep "HarmonicBondForce" log | awk '{ print $3 }' > pot_bond.log
	grep "HarmonicAngleForce" log | awk '{ print $3 }'  > pot_ang.log
	grep "PeriodicTorsionForce" log | awk '{ print $3 }' > pot_ptor.log
	grep "RBTorsionForce" log | awk '{ print $3 }' > pot_rtor.log
	grep "NonbondedForce" log | awk ' NR % 2 == 1' | awk '{ print $3 }' > pot_nbond.log
	grep "DrudeForce" log | awk '{ print $3 }' > pot_drude.log
	grep "CustomNonbondedForce" log | awk '{ print $3 }' > pot_cnbond.log
	grep "CustomBondForce" log | awk '{ print $3 }' > pot_cbond.log
fi
