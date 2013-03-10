#!/usr/bin/env python

__doc__ ="""
GREATx is an extension of the GREAT functionality. This module extends GREAT
by assigning weights to hits on the regulatory domains rather than counting
each hit equally. The classes and methods in this module use this model of
assigning weights to hits on the regulatory domain to calculate statistics.

Weights are determined for each dart in the regulatory domain by overlapping
normal distributions centered at each transcription start site. The regulatory
domain is defined as all bases within the defined cut-off distance of a single
transcription start site. Each instance of :WeightedRegDom: has the following
attributes:

    - sorted array of transcription start sites
    - cut-off distance to define the regulatory domain
    - mean value for the weight function (usually 0)
    - standard deviation for the weight function

With these attributes, each instance of :WeightedRegDom: can compute the weight
assigned to any dart on the genome and determine the maximum weight
assignment in the regulatory domain.

The different statistics included in this module are a Beta distribution
p-value, and the Gi-local, G general global, and Moran's I spatial
autocorrelation statistics.
"""

from scipy.stats import norm
from math import fabs
import collections
import re
import sys
import os
from collections import namedtuple


Result = namedtuple('Result', ['maxWgt', 'maxWgtDart'])

class WeightedRegDom:
    """
    Object which computes dart weights in a regulatory domain

    Parameters
    ----------
    transcriptionStartSites : array of transcription
                              start sites
    cutOff : maximum distance an dart can be from
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
    
    >>> twoWgt = wgtRegDom.getDartWgt(2)
    
    Get maximum site weight
    
    >>> maxWgt = wgtRegDom.maxDartWgt()
    
    """
    
    def __init__(self, transcriptionStartSites, cutOff, mean, sd):
        self.transcriptionStartSites = transcriptionStartSites
        self.transcriptionStartSites.sort()
        self.cutOff = cutOff
        self.wgtDistribution = self.initWgtDistribution(mean, sd, cutOff)

    def getWgtDistribution(self):
        """Return array of weight values for each base inside the cutOff for a tss."""
        return self.wgtDistribution

    def getTranscriptionStartSites(self):
        """Return tss array."""
        return self.transcriptionStartSites

    # note that this algorithm is okay because transcriptionStartSites is sorted
    def maxDartWgt(self):
        """Return namedtuple containing the max total weight and where it occurs."""

        # only do this hard calculation once!
        if (self.maxWgt):
            return Result(self.maxWgt, self.maxWgtDart)

        maxWgt = 0
        dart = 0
        for tss in self.transcriptionStartSites:
            for i in range(-self.cutOff, self.cutOff+1):
                if (dart < tss + i):  # dart only increases
                    dart = tss + i
                    dartWgt = self.getDartWgt(dart)
                    if (dartWgt > maxWgt):
                        maxWgt = dartWgt
                        maxWgtDart = dart
    
        return Result(self.maxWgt, self.maxWgtDart)

    def getDartWgt(self, dart):
        """Return the total weight assigned to an dart."""
        wgt = 0
        for tss in self.transcriptionStartSites:
            if (fabs(dart - tss) <= self.cutOff):
                wgt += self.getRelativeDartWgt(dart, tss)

        return wgt

    def getAllNonZeroDarts(self):
        """Return all darts which get a weight from a tss."""
        lastSite = self.transcriptionStartSites[-1]
        
        darts = []
        for tss in self.transcriptionStartSites:
            darts += [tss + i for i in range(-self.cutOff, self.cutOff+1)]
    
        return darts

    def getRelativeDartWgt(self, dart, tss):
        """Return the weight for a particular dart relative to a tss."""
        return self.wgtDistribution[tss-dart]

    def initWgtDistribution(self, mean, sd, cutOff):
        """Initialize array of weight values for each base inside the cutOff for a tss."""
        wgt = norm(mean, sd)
        maxWgt = wgt.pdf(mean)
        return [wgt.pdf(i)/maxWgt for i in range(-cutOff, cutOff + 1)]

    def readTSSPositions(fstr):
        """Read a file of numbers into an array"""
        tss = []
        f = open(fstr)
        for line in f:
            tss.append(int(line))
    
        return tss

def parseLociLine(lociLine, regSz):
    """Expands fields of lociLine to a regulatory region line."""
    lineFields = lociLine.split()
    chrName = lineFields[1]
    TSSPosition = lineFields[2]
    regStart = str(max(0, int(TSSPosition) - regSz))
    regEnd = str(int(TSSPosition) + regSz)
    geneName = lineFields[4]
    geneID = lineFields[0]
    strand = lineFields[3]

    return "\t".join([chrName, regStart, regEnd, geneName, geneID, strand, TSSPosition])
    #return line[1] +"\t"+str(max(0, int(line[2]) - regSz))+"\t"+str(int(line[2]) + regSz)+"\t"+line[4]+"\t"+line[0]+"\t"+line[3]+"\t"+line[2]+"\n")

def expandTSSsToRegDoms(lociFn, regDomFn, regSz):
    """Expands transcription start sites into regulatory regions"""
    loci = open(lociFn, 'r')
    regDom = open(regDomFn, 'w')
    #regSz = 1000000
    
    for line in loci:
        regDom.write(parseLociLine(line, regSz)+"\n")

def overlapSelectMergeOutput(regDomFn, dartFn, mergedFn):
    """Runs overlapselect -mergeoutput <regDomFn> <dartFn> <mergedFn>"""
    program = "/afs/ir/class/cs173/bin/i386_linux26/overlapSelect -mergeOutput"
    os.system(" ".join([program, regDomFn, dartFn, mergedFn]))

def getTSSs(lociFn):
    """Returns transcription start sites in a dictionary for each chromosome"""
    loci = open(lociFn, 'r')
    TSSs = {}

    for line in loci:
        line = line.split()
        chrName = line[1]
        TSSPosition = line[2]
        
        if (not (chrName in TSSs)):
            TSSs[chrName] = [TSSPosition]
        else:
            TSSs[chrName].append(TSSPosition)

    return TSSs
        
def assignWeights(cutOff, mean, sd, mergedFn, dartsToWeightsFn):
    """Writes to dartsToWeightsFn each dart with the geneName, geneID, and weight

     merged file must follow this format:
    |<-----------Dart Information---------->|<---------RegDom Info------->| Gene Name | GeneID | strand | TSS Location |
    | chrName | chrStart | chrEnd | dartName | chrName | chrStart | chrEnd | Gene Name | GeneID | strand | TSS Location |

    dartsToWeights file will be in this format:
    | dartName | Gene Name | Gene ID | Weight [0-1] |

    """

    merged = open(mergedFn, 'r')
    dartsToWeights = open(dartsToWeightsFn, 'w')

    for line in merged:
        line = line.split()
        dartPos = (int(line[1]) + int(line[2]))/2
        lociCenter = int(line[10])
        dartName = line[3] 
        geneName = line[7]
        geneID = line[8]

        wgtDist = norm(mean, sd)
        maxWgtDist = wgtDist.pdf(mean)
        dartWgt = wgtDist.pdf(lociCenter - dartPos)/maxWgtDist

        dartsToWeights.write("\t".join([dartName, geneName, geneID, repr(dartWgt)]) + "\n")

class Dart:

    def __init__(self, name, gene_name, id, weight):
        self.name = name
        self.gene_name = gene_name
        self.gene_id = id
        self.weight = weight

    def __str__(self):
        return 'Dart Name: ' + self.name + '\nGene Name: ' + self.gene_name + '\nID: ' + str(self.id) + '\nWeight: ' + str(self.weight)


class RegDom:
    """Instantiates a regulatory domain"""
    def __init__(self, start, end, id):
        self.id = id
        self.start = start
        self.end = end


class AssociationMaker:
    darts = []
    genetoterms = collections.defaultdict(lambda :[])
    termtocoverage = collections.defaultdict(lambda : 0.0)

    def __init__(self, dartsToWeightsFn, geneOntologyFn, regDomFn, chromSizesFn):
        self.readDartWeightsFile(dartsToWeightsFn)
        self.buildGeneTermMap(geneOntologyFn)
        self.buildTermWeightsMap(regDomFn, chromSizesFn)

    def readDartWeightsFile(self, fstr):
        with open(fstr) as f:
            for line in f:
                line = line.split("\t")
                dart = Dart(line[0], line[1], int(line[2]), float(line[3]))
                self.darts.append(dart)

    def buildGeneTermMap(self, geneOntologyFn):
        f = open(geneOntologyFn)
        for line in f:
            line = line.split("\t")
            term_id = int(re.search("\d+", line[0]).group(0))
            gene_id = int(re.search("\d+", line[1]).group(0))
            self.genetoterms[gene_id].append(term_id)

    def buildTermWeightsMap(self, regDomFn, chromSizesFn):
        genes = [dart.gene_id for dart in self.darts]
        regdoms = []
        with open(regDomFn) as f:
            for line in f:
                line = line.split()
                if int(line[4]) in genes:
                    regdoms.append(RegDom(int(line[1]),int(line[2]),int(line[4])))
        genome_size = 0
        with open(chromSizesFn) as f:
            for line in f:
                line = line.split()
                chrom_size = int(line[1])
                genome_size += chrom_size
        # build the map of terms to regdoms
        termtoregdoms = collections.defaultdict(lambda : [])
        for regdom in regdoms:
            terms = self.getTerms(regdom.id)
            for term in terms:
                termtoregdoms[term].append(regdom)
        # find the total length of all regions associated with that term
        for item in termtoregdoms.items():
            term = item[0]
            regdoms = item[1]
            nonoverlapping = self.removeOverlaps([(regdom.start,regdom.end) for regdom in regdoms])
            coverage = sum([x[1]-x[0] for x in nonoverlapping])
            percent = float(coverage)/genome_size
            self.termtocoverage[term] = percent


    # from http://stackoverflow.com/questions/1233292/whats-a-good-generic-algorithm-for-collapsing-a-set-of-potentially-overlapping
    def removeOverlaps(self, ranges):
        result = []
        cur = None
        for start, stop in sorted(ranges): # sorts by start
            if cur is None:
                cur = (start, stop)
                continue
            cStart, cStop = cur
            if start <= cStop:
                cur = (cStart, max(stop, cStop))
            else:
                result.append(cur)
                cur = (start, stop)
        result.append(cur)
        return result

    def getTerms(self, gene_id):
        #return [str(x[0]) for x in self.termtogenes.items() if gene in x[1]]
        return map(lambda x: str(x), self.genetoterms[gene_id])

    def buildLine(self, term_name, dart, term_weight):
        return "\t".join([term_name, str(dart.gene_id), dart.name, str(dart.weight), str(term_weight)]) + "\n"

    def writeOutput(self, output_file):
        f = open(output_file, "w")
        for dart in self.darts:
            terms = self.getTerms(dart.gene_id)
            if terms == []: # in case we get a gene that for some reason has no terms associated
                f.write(self.buildLine("UNKNOWN", dart, "0.0"))
            else:
                for term in terms:
                    f.write(self.buildLine(term, dart, self.termtocoverage[term]))
        f.close()

## DEPRECATED
if __name__ == '__main__':
    import sys

    from optparse import OptionParser
    parser = OptionParser(usage="%prog <tssFn> <cut-off> <mean> <standard deviation>",
                          description=("Determines the maximum possible weight "
                          "for darts landing in the regulatory domain "
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
    r = wgtRegDom.maxDartWgt()

    print repr(r.maxWgt) + "\t" + repr(r.maxWgtDart)

