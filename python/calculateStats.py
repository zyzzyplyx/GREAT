#!/usr/bin/env python

import sys
import os
import GREATx 
import scipy.stats
import collections

sys.path.append("./python")

## get data/SRFtoTerms.data
#GREATx.createRegDomsFileFromTSSs("GREATRegDoms/ontologies/hg18/hg18.loci", "data/hg18.regDom.bed", 1000000)
#GREATx.overlapSelect("data/hg18.regDom.bed", "GREATRegDoms/SRF.hg18.bed", "data/regDom.SRF.merge", options="-mergeOutput")
#GREATx.assignWeights(1000000, 0, 333333, "data/regDom.SRF.merge", "data/SRF.wgt")
#maker = GREATx.AssociationMaker("data/SRF.wgt", "GREATRegDoms/ontologies/hg18/GOBiologicalProcess/ontoToGene.canon", "data/hg18.regDom.bed", "GREATRegDoms/ontologies/hg18/chrom.sizes")
#maker.writeOutput("data/SRFtoTerms.data")

inFn = sys.argv[1]
outFn = sys.argv[2]
whichBeta = int(sys.argv[3])

inFile = open(inFn)
outFile = open(outFn, 'w')

weightSum = 0 #Will store our alpha
X = 0 #Will store our input probability 

lineObjects = []
for line in inFile:
    lineObjects.append(GREATx.TermDartTSSTriple(line))

termIDs = list(set([lineObject.termID for lineObject in lineObjects]))

termPValues = {}

for termID in termIDs:
    if termID != 'UNKNOWN':
        termIDObjects = filter(lambda x: x.termID == termID, lineObjects)
        weights = [termIDObject.weight for termIDObject in termIDObjects]

        alpha = sum(weights)
        x = termIDObjects[0].percentCoverage
        
        if whichBeta == 1:
            beta = len(termIDObjects) - alpha
        elif whichBeta == 2:
            beta = len(termIDObjects) * max(weights) - alpha
        elif whichBeta == 3:
            wgtRegDom = GREATx.WeightedRegDom(cutoff=1000000, mean=0, sd=333333)
            termTSSs = [termIDObject.TSSPosition for termIDObject in termIDObjects]
            beta = len(termIDObjects) * (wgtRegDom.bestWeightedDart(termTSSs, chromosomes=GREATx.HUMAN_CHROMOSOMES)).weight

        outFile.write(termID + "\t" + str(scipy.stats.beta.cdf(x, alpha, beta)) + "\n")

