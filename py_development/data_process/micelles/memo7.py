#!/home/kjeong23/softwares/bin/python3.4

import numpy
#preline="mic#  31 N=  29 r=   1.6322 A=   0.0250 sv2abc=   5.7960   4.0812   3.4426 members: [100, 112, 192, 203, 214, 226, 280, 295, 307, 339, 347, 396, 406, 420, 435, 439, 474, 588, 600, 610, 620, 622, 641, 649, 683, 799, 819, 849, 897]"

#memline=preline[84:]
#memline=memline.replace("]","")
#memline=memline.replace("[","")
#memline=memline.replace(",","")
#target=memline.split()
#print(target) 

#morline='step# 0 numofmcls 34\n'
#split=morline.split()
#print(split)
#print(split[3])

fcrd=numpy.array([[0,1,0],[4,0,0],[0,4,0]])
apex=numpy.array([1,8,6])

r12,r23,r31=fcrd[0]-fcrd[1],fcrd[1]-fcrd[2],fcrd[2]-fcrd[0]
sa,sb,sc=numpy.linalg.norm(r12),numpy.linalg.norm(r23),numpy.linalg.norm(r31)
ss=(sa+sb+sc)/2.0
oneS=numpy.sqrt( ss * (ss-sa) * (ss-sb) * (ss-sc) )

#Volume of trigonal pyramid is (1/3)*A*H.
pnorm=numpy.cross(r23,r12)
#pnorm=pnorm/( numpy.linalg.norm(pnorm) )
r14=fcrd[0]-apex
h=numpy.linalg.norm( numpy.dot(pnorm,r14) ) /numpy.linalg.norm(pnorm)
oneV=oneS*h/3.0

SandV=numpy.array([oneS,oneV])
print (pnorm)
print (SandV)

