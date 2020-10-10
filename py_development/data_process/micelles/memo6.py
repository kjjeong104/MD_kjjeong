#!/home/kjeong23/softwares/bin/python3.4

import numpy

#hopt=numpy.array([[150,300,550],[90,-1,-1],[50,-1,-1],[234,500,-1]])

#print(hopt)
#print(numpy.argwhere(hopt[1]==-1))

#updpos=numpy.argwhere(hopt[1]==-1)
#hopt[1][updpos[0]]=340
#print(hopt)

#nsurf=10
#unagg=[x for x in range(nsurf)]
#print(unagg)
#unagg.remove(7)
#print(unagg)
#unagg.remove(9)
#print(unagg)

empty=numpy.empty((3,0),float)
attach=[-1]*3
empty=numpy.hstack((empty,attach))
print(empty)
