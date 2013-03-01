#!/usr/bin/env python
#Usage: ./findDist.py overlapFile outputFile

import os, sys, scipy.stats

overlap = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')
regSz = 1000000

dist = scipy.stats.norm(0,regSz / 3)
w_max = dist.pdf(0)

for line in overlap:
    line = line.split()
    arrowPos = int(line[1])
    lociCenter = int(line[10])

    distance = abs(lociCenter - arrowPos)

    outFile.write(line[3] + "\t" +  line[7] + "\t" + line[8] + "\t" +
    str(dist.pdf(distance) / w_max) + "\n")

overlap.close()
outFile.close()
