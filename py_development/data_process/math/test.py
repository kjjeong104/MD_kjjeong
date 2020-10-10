#!/home/kjeong23/softwares/bin/python3.4

import math
import numpy
import numpy.linalg

xvec=numpy.array([1,2,3,4,5])
xvec=numpy.append(xvec,6)
xvec=numpy.delete(xvec,0)
print(xvec)

B=numpy.array([[1,2],[3,4]])
cvec=numpy.array([0,1])
print(B.dot(cvec))
