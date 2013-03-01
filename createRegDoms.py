#!/usr/bin/env python
#Usage: ./createRegDoms.py lociFile outputFile

import os, sys

loci = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')
regSz = 1000000

for line in loci:
  line = line.split()
  #print line[1]+" 1: "+line[1]
  #TODO: Change this to cap the region at the end of the chromosome
  outFile.write(line[1] +"\t"+str(max(0, int(line[2]) - regSz))+"\t"+str(int(line[2]) + regSz)+"\t"+line[3]+"\t"+line[4]+"\t"+line[0]+"\t"+line[2]+"\n")


  
