#!/bin/sh
num=1
#mkdir series_test

while [ $num -le 100 ]
do
	./water_spc
	cp spc_water.dat ./series_test/spc_water$num.dat
	cp old_config.ini ./series_test/old_config$num.ini
	cp config_pov.dat ./series_test/config_pov$num.dat
	num=`expr $num + 1`
done

