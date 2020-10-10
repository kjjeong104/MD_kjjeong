import numpy as np
import sys

set1 = np.loadtxt(sys.argv[1]+'/pot_sum_nbond.log')
set2 = np.loadtxt(sys.argv[2]+'/pot_sum_nbond.log')
data = set1 - set2
np.savetxt('pot_sum_nbond_intra.log',data)
np.savetxt('pot_sum_nbond_inter.log',set2)

