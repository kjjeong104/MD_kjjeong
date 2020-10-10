#!/home/kjeong23/softwares/bin/python3.4

import numpy

array1=[1,2,3, \
4,5,6]
array2=numpy.array([5,6,7,8,9,10])
array5=[9]

array3=[[3,4],[1,2],[5,6],[7,8]]
array4=[[9,10],[11,12],[13,14],[15,16]]

#print (set(array1+array2+array5))
#print (array3[0]+array4[3])

array3.sort(reverse=True)
array3.append([2,1,4,array2[3]])
del array3[0]
print(array3)

print('numpy test')
print(array2)
array22 = array2 / 2
print(array22)
#array22=numpy.append(array2,4)
print(array22-array2)
array6=numpy.append(array2,[1,2,3,'one','two','three'])
print(array6)
print(array6.size)
#for i in array3:
#  print (i)

#if [11,12] in array4:
#  print('1d array [11,12] is detected in array4')

#if [17,18] not in array4:
#  print('[17,18] is not in array4')

#for i in range(len(array2)):
#  for j in range(i+1,len(array2)):
#    print (i,j)
