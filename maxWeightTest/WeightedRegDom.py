#!/usr/bin/env python

__doc__ ="""
Weights are determined for each arrow in the regulatory domain by overlapping
normal distributions centered at each transcription start site. The regulatory
domain is defined as all bases within the defined cut-off distance of a single
transcription start site. Each instance of :WeightedRegDom: has the following
attributes:
    - sorted array of transcription start sites
    - cut-off distance to define the regulatory domain
    - mean value for the weight function (usually 0)
    - standard deviation for the weight function

With these attributes each instance of :WeightedRegDom: can compute the weight
assigned to any arrow on the genome and determine the maximum weight
assignment in the regulatory domain.
"""

from scipy.stats import norm as _norm
from math import fabs as _fabs
from collections import namedtuple

_Result = namedtuple('Result', ['maxWgt', 'maxWgtArrow'])

class WeightedRegDom:

    """
    Object which computes arrow weights in a regulatory domain

    Parameters
    ----------
    transcriptionStartSites : array of transcription
                              start sites
    cutOff : maximum distance an arrow can be from
             a transcription start site to be
             assigned any weight from it
    mean : mean value of the normal distribution
           for calculating weights
    sd : standard deviation of the normal distribution
         for calculating weights

    
    Examples
    --------
    >>> import WeightedRegDom as wrd
    >>> tss = [1, 2, 3]
    >>> wgtRegDom = wrd.WeightedRegDom(tss, cutOff=2, mean=0.0, sd=3)
    
    Get weight for a given site
    
    >>> twoWgt = wgtRegDom.getArrowWgt(2)
    
    Get maximum site weight
    
    >>> maxWgt = wgtRegDom.maxArrowWgt()
    
    """
    
    def __init__(self, transcriptionStartSites, cutOff, mean, sd):
        self.transcriptionStartSites = transcriptionStartSites
        self.transcriptionStartSites.sort()
        self.cutOff = cutOff
        self.wgtDistribution = self._initWgtDistribution(mean, sd, cutOff)

    def getWgtDistribution(self):
        """Return array of weight values for each base inside the cutOff for a tss."""
        return self.wgtDistribution

    def getTranscriptionStartSites(self):
        """Return tss array."""
        return self.transcriptionStartSites

    # note that this algorithm is okay because transcriptionStartSites is sorted
    def maxArrowWgt(self):
        """Return namedtuple containing the max total weight and where it occurs."""
        maxWgt = 0
        arrow = 0
        for tss in self.transcriptionStartSites:
            for i in range(-self.cutOff, self.cutOff+1):
                if (arrow < tss + i):  # arrow only increases
                    arrow = tss + i
                    arrowWgt = self.getArrowWgt(arrow)
                    if (arrowWgt > maxWgt):
                        maxWgt = arrowWgt
                        maxWgtArrow = arrow
    
        return _Result(maxWgt, maxWgtArrow)

    def getArrowWgt(self, arrow):
        """Return the total weight assigned to an arrow."""
        wgt = 0
        for tss in self.transcriptionStartSites:
            if (_fabs(arrow - tss) <= self.cutOff):
                wgt += self._getRelativeArrowWgt(arrow, tss)

        return wgt

    def getAllNonZeroArrows(self):
        """Return all arrows which get a weight from a tss."""
        lastSite = self.transcriptionStartSites[-1]
        
        arrows = []
        for tss in self.transcriptionStartSites:
            arrows += [tss + i for i in range(-self.cutOff, self.cutOff+1)]
    
        return arrows

    def _getRelativeArrowWgt(self, arrow, tss):
        """Return the weight for a particular arrow relative to a tss."""
        return self.wgtDistribution[self.cutOff-(tss-arrow)]


    def _initWgtDistribution(self, mean, sd, cutOff):
        """Initialize array of weight values for each base inside the cutOff for a tss."""
        wgt = _norm(mean, sd)
        maxWgt = wgt.pdf(mean)
        return [wgt.pdf(i)/maxWgt for i in range(-cutOff, cutOff + 1)]

def readPositions(fstr):
    """Read a file of numbers into an array"""
    tss = []
    f = open(fstr)
    for line in f:
        tss.append(int(line))

    return tss

if __name__ == '__main__':
    import sys

    from optparse import OptionParser
    parser = OptionParser(usage="%prog <tssFn> <cut-off> <mean> <standard deviation>",
                          description=("Determines the maximum possible weight "
                          "for arrows landing in the regulatory domain "
                          "defined by <tssFn> and <cut-off> with "
                          "overlapping normal distributions shaped by <mean> and "
                          "<standard deviation>."
                          "The file <tssFn> should have a number corresponding "
                          "to each transcription start site on every line."))
    #parser.add_option("-q", "--quiet",
    #                  action="store_false", dest="verbose", default=True,
    #                  help="don't print status messages to stdout")
    
    (options, args) = parser.parse_args()
    if (len(args) != 4):
        parser.print_usage()
        sys.exit(1)

    tss    = readPositions(args[0])
    cutOff = int(args[1])
    mean   = float(args[2])
    sd     = float(args[3])

    wgtRegDom = WeightedRegDom(tss, cutOff, mean, sd)
    r = wgtRegDom.maxArrowWgt()

    print repr(r.maxWgt) + "\t" + repr(r.maxWgtArrow)

