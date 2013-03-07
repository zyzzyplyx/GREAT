#!/usr/bin/env python

#Usage: ./maxWeight.py transcriptionStartSites

import sys, scipy.stats, math, collections

Result = collections.namedtuple('Result', ['maxWgt', 'maxWgtMarker'])

class NormWgt:

    def __init__(self, mean, sd, cutOff, transcriptionStartSites):
        self.transcriptionStartSites = transcriptionStartSites
        self.cutOff = cutOff
        self.wgtDistribution = self.getWgtDistribution(mean, sd, cutOff)
        print repr(self.wgtDistribution) + "\n\n"

    def getSiteWgtTotal(self, site):
        wgtTotal = 0
        for tss in self.transcriptionStartSites:
            if (math.fabs(site - tss) <= self.cutOff):
                wgtTotal += self.getSiteWgt(site, tss)

        return wgtTotal

    def getSiteWgt(self, site, tss):
        return self.wgtDistribution[self.cutOff-(tss-site+1)]

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
    normWgt = NormWgt(MEAN, SD, CUTOFF, transcriptionStartSites)
    nonZeroMarkers = getAllNonZeroMarkers(transcriptionStartSites, CUTOFF)

    maxWgt = 0
    for marker in nonZeroMarkers:
        markerWgt = normWgt.getSiteWgtTotal(marker)
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
    MEAN = float(sys.argv[2])
    SD = float(sys.argv[3])
    CUTOFF = int(sys.argv[4])

    r = findMaxNormPDFOverlap(MEAN, SD, CUTOFF, transcriptionStartSites)
    print repr(r.maxWgt) + "\t" + repr(r.maxWgtMarker)

if __name__ == '__main__':
    main()
