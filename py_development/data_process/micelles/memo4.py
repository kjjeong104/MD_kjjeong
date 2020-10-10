#!/home/kjeong23/softwares/bin/python3.4

import numpy
import numpy.linalg

#a1=numpy.array([5,6,7,8,9,10])
#print (a1.size)
a2=numpy.array([[1,5],[2,3],[1,6]])
#a3=numpy.array([])
print (a2)
print (a2[0:2])
print (a2.shape[0])
print (a2.shape[1])
#for i in a2:
#  print(i)
#a4=numpy.array([5,12])
#print (numpy.sqrt(numpy.dot(a4,a4)))
#a2=numpy.zeros((3,3))
#print(a2)
#print(a2[2,2])
#print(a2[0])
#a3=numpy.array([[1,2,3],[4,5,6],[7,8,9]])
#print(a3)
#e_values, e_vectors = numpy.linalg.eig(a3)
#print(e_values)
a5=numpy.array([[1,1],[1,1],[1,1]])
print(a2+a5)

print(numpy.mean(a2,axis=0))
