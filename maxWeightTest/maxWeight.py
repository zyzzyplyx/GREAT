#!/usr/bin/env python

#Usage: ./maxWeight.py transcriptionStartSites MEAN SD CUTOFF

import sys, scipy.stats, math, collections

Result = collections.namedtuple('Result', ['maxWgt', 'maxWgtMarker'])

class NormWgt:

    def __init__(self, mean, sd, cutOff, transcriptionStartSites):
        self.transcriptionStartSites = transcriptionStartSites
        self.cutOff = cutOff
        self.wgtDistribution = self.getWgtDistribution(mean, sd, cutOff)

    def getSiteWgtTotal(self, site):
        wgtTotal = 0
        for tss in self.transcriptionStartSites:
            if (math.fabs(site - tss) <= self.cutOff):
                wgtTotal += self.getSiteWgt(site, tss)

        return wgtTotal

    def getSiteWgt(self, site, tss):
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

### this is an idea that is poorly formed. it feels like a waste to look at
## every base when looking for the max
#
#def getAllCandidateMarkers(transcriptionStartSites, CUTOFF):
#    transcriptionStartSites.sort()
#    lastSite = transcriptionStartSites[-1]
#    firstSite = transcriptionStartSites[0]
#    lastIndex = len(transcriptionStartSites) - 1
#    
#    candidateMarkers = []
#    beforeIndex = 0
#    for i in range(max(0,firstSite-CUTOFF), lastSite+CUTOFF+1):
#        if (beforeIndex <= lastIndex and i == transcriptionStartSites[beforeIndex]):
#            beforeIndex += 1
#            candidateMarkers.append(i)
#        else:
#            if (beforeIndex == 0):
#                c = int(transcriptionStartSites[beforeIndex] - i <= CUTOFF)
#                if (beforeIndex + 1 <= lastIndex):
#                    d = int(transcriptionStartSites[beforeIndex+1] - i <= CUTOFF)
#                else:
#                    d = 0
#                if (c+d >= 2):
#                    candidateMarkers.append(i)
#            elif (beforeIndex == 1):
#                b = int(i - transcriptionStartSites[beforeIndex-1] <= CUTOFF)
#                c = int(transcriptionStartSites[beforeIndex] - i <= CUTOFF)
#                if (beforeIndex + 1 <= lastIndex):
#                    d = int(transcriptionStartSites[beforeIndex+1] - i <= CUTOFF)
#                else:
#                    d = 0
#                if (b+c+d >= 2):
#                    candidateMarkers.append(i)
#            elif (beforeIndex == lastIndex):
#                if (beforeIndex - 2 >= 0):
#                    a = int(i - transcriptionStartSites[beforeIndex-2] <= CUTOFF)
#                else:
#                    a = 0
#                if (beforeIndex - 1 >= 0):
#                    b = int(i - transcriptionStartSites[beforeIndex-1] <= CUTOFF)
#                else:
#                    b = 0
#                c = int(transcriptionStartSites[beforeIndex] - i <= CUTOFF)
#                if (a+b+c >= 2):
#                    candidateMarkers.append(i)
#            elif (beforeIndex > lastIndex):
#                if (beforeIndex - 2 >= 0):
#                    a = int(transcriptionStartSites[beforeIndex-2] - i <= CUTOFF)
#                else:
#                    a = 0
#                b = int(transcriptionStartSites[beforeIndex-1] - i <= CUTOFF)
#                if (a+b >= 2):
#                    candidateMarkers.append(i)
#            else:
#                a = int(i - transcriptionStartSites[beforeIndex-2] <= CUTOFF)
#                b = int(i - transcriptionStartSites[beforeIndex-1] <= CUTOFF)
#                c = int(transcriptionStartSites[beforeIndex] - i <= CUTOFF)
#                d = int(transcriptionStartSites[beforeIndex+1] - i <= CUTOFF)
#                if (a+b+c+d >= 2):
#                    candidateMarkers.append(i)
#
#    return candidateMarkers

def findMaxNormPDFOverlap(MEAN, SD, CUTOFF, transcriptionStartSites):
    normWgt = NormWgt(MEAN, SD, CUTOFF, transcriptionStartSites)
    nonZeroMarkers = getAllNonZeroMarkers(transcriptionStartSites, CUTOFF)
    #nonZeroMarkers = getAllCandidateMarkers(transcriptionStartSites, CUTOFF)

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
