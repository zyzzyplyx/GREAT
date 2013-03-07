#!/usr/bin/env python

#Usage: ./maxWeight.py transcriptionStartSites

import sys, scipy.stats, math, collections

MEAN = 0
SD = 10**1/3
CUTOFF = 10**1

Result = collections.namedtuple('Result', ['maxWgt', 'maxWgtMarker'])

class NormPull:

    def __init__(self, mean, sd, cutOff, transcriptionStartSites):
        self.transcriptionStartSites = transcriptionStartSites
        self.cutOff = cutOff
        self.wgtDistribution = self.getWgtDistribution(mean, sd, cutOff)

    def getSitePullTotal(self, site):
        pullTotal = 0
        for tss in self.transcriptionStartSites:
            if (math.fabs(site - tss) <= self.cutOff):
                pullTotal += self.getSitePull(site, tss)

        return pullTotal

    def getSitePull(self, site, tss):
        return self.wgtDistribution[self.cutOff-(tss-site)]

    def getWgtDistribution(self, mean, sd, cutOff):
        wgt = scipy.stats.norm(mean, sd)
        return [wgt.pdf(i) for i in range(-cutOff, cutOff + 1)]

def getAllNonZeroMarkers(transcriptionStartSites, CUTOFF):
    transcriptionStartSites.sort()
    lastSite = transcriptionStartSites[-1]
    
    markers = []
    for tss in transcriptionStartSites:
        markers += [tss + i for i in range(-CUTOFF, CUTOFF+1)]

    return markers

def findMaxNormPDFOverlap(MEAN, SD, CUTOFF, transcriptionStartSites):
    normPull = NormPull(MEAN, SD, CUTOFF, transcriptionStartSites)
    nonZeroMarkers = getAllNonZeroMarkers(transcriptionStartSites, CUTOFF)

    maxWgt = 0
    for marker in nonZeroMarkers:
        markerWgt = normPull.getSitePullTotal(marker)
        if (markerWgt > maxWgt):
            maxWgt = markerWgt
            maxWgtMarker = marker

    return Result(maxWgt, maxWgtMarker)

def readPositions(fstr):
    tss = []
    f = open(fstr)
    for line in f:
        tss.append(int(line))
        
    return tss

def main():
    transcriptionStartSites = readPositions(sys.argv[1])
    r = findMaxNormPDFOverlap(MEAN, SD, CUTOFF, transcriptionStartSites)
    print repr(r.maxWgt) + "\t" + repr(r.maxWgtMarker)

if __name__ == '__main__':
    main()
