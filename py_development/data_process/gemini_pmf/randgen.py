#!/home/kjeong23/softwares/bin/python3.4

import math
import sys
import numpy
import random
from numpy import linalg

def main():

  for i in range(22):
    print("dummy line")

  for i in range(1000):
    x,y,z=(random.random()-0.5)*2.0,(random.random()-0.5)*2.0,(random.random()-0.5)*2.0
    r=numpy.linalg.norm(numpy.array([x,y,z]))
    print("{:8.4f} {:8.4f} {:8.4f} {:8.4f} {:8.4f}".format(i*1.0,r,x,y,z))

if __name__ == "__main__": main()


