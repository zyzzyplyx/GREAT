#!/usr/bin/env python
#
# Usage: ./findDist.py overlapFile outputFile
#
# Takes a list of overlapping "arrows" and their corresponding region,
# calculates the weight for each overlap. Uses a normal distribution centered
# at the location of the TSS with std dev of 333,333 bp.
#
# overlapFile must follow this format:
#|<--------------Arrow Information--------------->|<-------RegDom Info---------->| Gene Name |  GeneID | strand | TSS Location |
#| chr5    60663738        60663788        SRF.1  | chr5   59031713    61031713  |  DEPDC1B  |  19474  |    -   |  60031713    |
#
# outputFile will be in this format:
# | “Arrow” name (e.g. SRF.1) | Gene Name | Gene ID (numerical ID) | Weight (in the range of 0 to 1) |
# |           SRF.1           |  DEPDC1B  |        19474           |          0.165703737676         |

import os, sys, scipy.stats

overlap = open(sys.argv[1], 'r')
outFile = open(sys.argv[2], 'w')
regSz = 1000000

dist = scipy.stats.norm(0,regSz / 3)
w_max = dist.pdf(0)

for line in overlap:
  line = line.split()
  arrowPos = (int(line[1]) + int(line[2]))/2
  lociCenter = int(line[10])

  distance = abs(lociCenter - arrowPos)

  outFile.write(line[3] + "\t" +  line[7] + "\t" + line[8] + "\t" +
  str(dist.pdf(distance) / w_max) + "\n")

overlap.close()
outFile.close()