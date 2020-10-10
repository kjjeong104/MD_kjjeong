#!/home/kjeong23/softwares/bin/python3.4

#format example for .gro:
#    1DEP     C1    1   1.748  15.989   1.574 -0.6456  0.9631 -0.0749
test=[1, 'DEP', 'X', 1, 1.748, 15.989, 1.574]
print ('{:5}{:5}{:5}{:5}{:8}{:8}{:8}'.format(test[0],test[1],test[2],test[3],test[4],test[5],test[6]))

test2='    1DEP      X    1   1.748  15.989   1.574 -0.6456  0.9631 -0.0749'
splitstr=[test2[0:10],test2[10:15],test2[15:20],test2[20:28],test2[29:36],test2[37:44]]
for i in range(len(splitstr)):
  splitstr[i]=splitstr[i].replace(" ","")
print(splitstr)
#str1=str(1)
#str2='DEP'
#str3=str1+str2
#print (str1+str2)
