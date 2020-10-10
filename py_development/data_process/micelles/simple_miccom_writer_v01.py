#!/home/kjeong23/softwares/bin/python3.4
# informal micelle com position with site(micelle) tracking.
# input file : merged micelle COM .gro file, 
#and tracker matrix file(made by kinmic_v04.py)
# simply read tracker matrix file, and according to reassigned site number,
# (considered position relevance, micelle movement history tracked in tracker matrix)
# write some informal file notates COM positions of well-assigned micelles

import sys
import numpy
import timeit

def com_rewrite(comstr,num,t): #rewrite COM position line with correcting the site number
  comsplit=comstr.split()
  #newstr="nlatp# "+'{:2}'.format(num)+str[9:]
  newstr="nlatp# "+'{:2} {:8.3f} ns '.format(num,t)+\
  '{:8} {:8} {:8}\n'.format(comsplit[3],comsplit[4],comsplit[5])
  return newstr

def main():
  #Load input files
  comfile = open(sys.argv[1],'r')
  matrixfile = open(sys.argv[2],'r')
  outfile=open(sys.argv[3],'w') #output file name. gave up to split files.

  #nnmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  start_time=timeit.default_timer()

  sindex=0
  while True:
    comline=comfile.readline()
    if comline[0:7]=="Centers":#step initiating line, start inner loop of 1step process
      paragraph=[]
      ltsplit1=comline.split()
      t=float(ltsplit1[6])
      comline=comfile.readline()
      nmcl=int(comline)
      #read 1step morphology
      for i in range(nmcl):
        comstr=comfile.readline()
        paragraph.append(comstr)
      dummy=comfile.readline() #boxlength line
      #read tracker matrix
      matline=matrixfile.readline()
      matsplit=matline.split()
      undismcl=nmcl
      for n in matsplit: #unassignable micelle test
        if n=='X':
          undismcl-=1
      j=0
      outfile.write('step# {} undislocated_nmcl {}\n'.format(sindex,undismcl))
      for n in matsplit:
        if n!='X':
          newline=com_rewrite(paragraph[j],n,t)
          outfile.write(newline)
        j+=1

      sindex+=1 #at the end of 1step loop, add step number
    if (sindex%500)==0:                               
      elapsed=timeit.default_timer() - start_time
      print('rewrote step {}, time {:11.4f}'.format(sindex,elapsed))
    if comline=='': #EOF
      break

  comfile.close()
  matrixfile.close()
  outfile.close()
   
if __name__ == "__main__": main()
