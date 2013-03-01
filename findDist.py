#!/usr/bin/env python
#Usage: ./findDist.py overlapFile outputFile

import os, sys

overlap = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')
regSz = 1000000


for line in overlap:
    line = line.split()
    arrowPos = int(line[1])
    lociCenter = int(line[10])

    distance = abs(lociCenter - arrowPos)

    outFile.write(line[3] + "\t" +  line[7] + "\t" + line[8] + "\t" +
    str(distance) + "\n")

overlap.close()
outFile.close()
