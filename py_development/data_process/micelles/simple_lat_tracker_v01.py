#!/home/kjeong23/softwares/bin/python3.4
# simple lattice site track-rewriter
# input file : merged lattice point assignment summary file, 
#and tracker matrix file(made by kinmic_v04.py)
# simply read tracker matrix file, and according to reassigned site number,
# (considered position relevance, micelle movement history tracked in tracker matrix)
## split the lattice morphology information into individual files by site. -> abandoned idea

import sys
import numpy
import timeit

def morph_rewrite(str,num): #rewrite morphology line with correcting the site number
  newstr="nlatp# "+'{:2}'.format(num)+str[9:]
  return newstr

def main():
  #Load input files
  morfile = open(sys.argv[1],'r')
  matrixfile = open(sys.argv[2],'r')
  outfile=open(sys.argv[3],'w') #output file name. gave up to split files.

  #nnmic=int(input("Most frequent number of micelles desired? ex) 30 or 8\n"))
  start_time=timeit.default_timer()

  sindex=0
  while True:
    morline=morfile.readline()
    if morline[0:4]=="step":#step initiating line, start inner loop of 1step process
      paragraph=[]
      ltsplit=morline.split()
      nmcl=int(ltsplit[3])
      #read 1step morphology
      for i in range(nmcl):
        morstr=morfile.readline()
        paragraph.append(morstr)
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
          newline=morph_rewrite(paragraph[j],n)
          outfile.write(newline)
        j+=1

      sindex+=1 #at the end of 1step loop, add step number
    if (sindex%200)==0:                               
      elapsed=timeit.default_timer() - start_time
      print('rewrote step {}, time {:11.4f}'.format(sindex,elapsed))
    if morline=='': #EOF
      break

  morfile.close()
  matrixfile.close()
  outfile.close()
   
if __name__ == "__main__": main()
