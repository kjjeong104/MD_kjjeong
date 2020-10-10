# you can select whenever you want
import sys
import numpy as np

if sys.argv[1] == 'no':
	set1 = np.loadtxt('pot_ang.log')
	set2 = np.loadtxt('pot_ptor.log')
	set3 = np.loadtxt('pot_bond.log')
	set4 = np.loadtxt('pot_rtor.log')
	data = (set1+set2+set3+set4)
	np.savetxt('pot_sum_bond.log',data)

	set1 = np.loadtxt('pot_cnbond.log')
	set2 = np.loadtxt('pot_nbond.log')
	set3 = np.loadtxt('pot_cbond.log')
	set4 = np.loadtxt('pot_drude.log')
	data = (set1+set2+set3+set4)
	np.savetxt('pot_sum_nbond.log',data)
elif sys.argv[1] == 'coh':
	set1 = np.loadtxt('pot_cnbond.log')
	set2 = np.loadtxt('pot_nbond.log')
	data = (set1+set2)
	np.savetxt('pot_sum_nbond.log',data)
else:
	print("wrong argument (coh/no)")
	exit


