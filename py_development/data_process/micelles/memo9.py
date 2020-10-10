#!/home/kjeong23/softwares/bin/python3.4

#format example for .gro:
#    1DEP     C1    1   1.748  15.989   1.574 -0.6456  0.9631 -0.0749
list=[3,4,5,6]
#print(list)
str1=str(list)+'{:8.3f}'.format(list[1])
print(str1)

morline=str1.replace("["," ")
morline=morline.replace("]",",")
print(morline)
