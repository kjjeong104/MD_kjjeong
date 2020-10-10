#!/home/kjeong23/softwares/bin/python3.4

import numpy
#line1="time    3.700 ns to    3.720 ns. accepted 0 donated 1  accmem [] donmem ['838']"
line1="time  112.820 ns to  112.840 ns. accepted 0 donated 3  accmem [] donmem ['11', '188', '451']"
#line1="time   47.060 ns to   47.080 ns. accepted 1 donated 0  accmem ['2'] donmem []"
sub1=line1[line1.index('donmem')+len('donmem'):]
sub1=sub1.replace("[","")
sub1=sub1.replace("]","")
sub1=sub1.replace("'","")
sub1=sub1.replace(",","")
split1=sub1.split()
#split1=line1.split('donmem')
#print(split1)
#print(split1)
#print(len(split1))

donlist=numpy.empty((0,2),float)
i=1.5
for x in split1:
  numrow=numpy.array([i,int(x)])
  donlist=numpy.vstack((donlist,numrow))
  i+=2.0

print(donlist)
#print(donlist.shape[0])
time_vector=donlist[:,0]
member_vector=donlist[:,1]
print(len(member_vector))
#print(time_vector)
#print(member_vector)
t1=2.7
i1=numpy.searchsorted(time_vector,t1,side='right')
print(i1)

