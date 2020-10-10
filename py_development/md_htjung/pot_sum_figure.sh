# script to collect potential components
#  for paper figures
# $1 : del/no (delete fist half values/no)

delete=$1

for folder in "gas_longMD" "liq_longMD"; do
	echo " move to " $folder
	cd $folder
	~/des/pot_component.sh no
	python ~/des/pot_sum.py no
	cd ..
done

for folder in "gas_longMD_coh" "liq_longMD_coh"; do
        echo " move to " $folder
        cd $folder
        ~/des/pot_component.sh coh
	python ~/des/pot_sum.py coh
        cd ..
done

# copy files
cp gas_longMD/pot_sum_bond.log ./pot_sum_bond_gas.log
cp liq_longMD/pot_sum_bond.log ./pot_sum_bond_liq.log

# nbond intra, inter
for prefold in "gas" "liq"; do
	python ~/des/pot_sum_nbond.py $prefold"_longMD" $prefold"_longMD_coh"
	prefile1="pot_sum_nbond_intra"
	prefile2="pot_sum_nbond_inter"
	mv $prefile1".log" $prefile1"_"$prefold".log"
	mv $prefile2".log" $prefile2"_"$prefold".log"

done

# delete?
#if [ $delete == 'del' ]; then
#	nlines=$(wc *.log | head -1 | awk '{ print $1 }')
#	hnlines=$($nlines/2)
#	sed -i "1,${hnlines}d" *.log
#fi

nlines=$(wc *.log | head -1 | awk '{ print $1 }')
rlines=$( expr $nlines % 2 )
if [ $rlines == 1 ]; then
	sed -i '1d' *.log
fi

# block average
for ifile in ./*.log; do
	if [ $delete == 'del' ]; then
		python ~/Utility/python/simple/blockavg.py -i $ifile -b -1
	else
		python ~/Utility/python/simple/blockavg.py -i $ifile -b 0
	fi
done

