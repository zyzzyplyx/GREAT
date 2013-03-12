#! /usr/bin/python2.7
"""GREAT extension classes

GREATx extends GREAT by assigning weighted hits on regulatory domains on
the continuous interval [0,1]. By assigning weights in this way,
we are able to generalize GREAT's Binomial distribution p-value to a
Beta distribution p-value. In addition, we add the Getis-Ord Gi* local
statistic -- a spatial autocorrelation statistic.
"""

__author__ = """Ahmed Bou-Rabee <bourabee@stanford.edu>,
                Andrew Stiles <aostiles@stanford.edu>,
                Charles Celerier <cceleri@cs.stanford.edu>, 
                Nicholas Damien McGee <ndmcgee@cs.stanford.edu>
                Steven Lee <slee2010@stanford.edu>, and
                Yifei Men <ymen@stanford.edu>"""
            
__date__ = """13 March 2012"""

__credits__ = """Gill Bejerano, for a thrilling tour of the genome.
Jim Notwell and Harendra Guturu, for their advice on this project."""

from scipy.stats import norm
from math import fabs
import collections
import re
import sys
import os
import operator
import scipy.stats

# these are the hard-coded human chromosome names and sizes
HUMAN_CHROMOSOMES = ['chr' + str(i) for i in range(1,23)] + ['chrX', 'chrY']
HUMAN_CHROMOSOME_SIZES = [ 247249719,\
                           242951149,\
                           199501827,\
                           191273063,\
                           180857866,\
                           170899992,\
                           158821424,\
                           146274826,\
                           140273252,\
                           135374737,\
                           134452384,\
                           132349534,\
                           114142980,\
                           106368585,\
                           100338915,\
                           88827254,\
                           78774742,\
                           76117153,\
                           63811651,\
                           62435964,\
                           46944323,\
                           49691432,\
                           154913754,\
                           57772954]

class Dart:
    """Object representing a dart on the human genome
    
    Darts are used to represent the darts placed on the genome which GREATx
    uses for statistical analysis based on those darts' "closeness" to the
    genes for each term in the list of GREAT ontologies.

    Parameters
    ----------
    chrName : str
              Name of the chromosome where the dart lives. (default = '')
    name : str
           Unique dart name. (default = '')
    position : int
               Unique dart position. (default = -1)

    Attributes
    ----------
    Same as parameters.

    Example
    --------
    >>> from GREATx import Dart
    >>> Dart(chrName='chr1', name='SRF.1', position='1000')
    GREATx.Dart(chrName="'chr1'", name="'SRF.1'", position='1000')
    """

    def __init__(self, **kwargs):
        self.chrName = kwargs.get('chrName', '')
        self.name = kwargs.get('name', '')
        self.position = int(kwargs.get('position', -1))

    def __str__(self):
        return "Dart Name: %s\nName: %s\nPosition: %d" % \
                (self.chrName, self.name, self.position)

    def __repr__(self):
        return 'GREATx.Dart(chrName=%r, name=%r, position=%r)' % \
                (repr(self.chrName), repr(self.name), repr(self.position))

class WeightedDart(Dart):
    """Object representing a dart on the human genome with a weight.
    
    A WeightedDart is a Dart with an addition attribute of weight. This
    object is particularly useful when calculating the Gi local statistic.

    Parameters
    ----------
    chrName : str
              Name of the chromosome where the dart lives. (default = '')
    name : str
           Unique dart name. (default = '')
    position : int
               Unique dart position. (default = -1)
    weight : float
             Assigned weight. (default = -1)       

    Attributes
    ----------
    Same as parameters.

    Example
    --------
    >>> from GREATx import WeightedDart
    >>> WeightedDart(chrName='chr1', name='SRF.1', position=1000, weight=23)
    GREATx.WeightedDart(chrName="'chr1'", name="'SRF.1'", position='1000', weight='23')
    """

    def __init__(self, **kwargs):
        self.weight = kwargs.get('weight', -1)
        super(WeightedDart, self).__init__(**kwargs)

    def __str__(self):
        return Dart.__str__(self) + ("\nWeight: %f" % self.weight)

    def __repr__(self):
        return 'GREATx.WeightedDart(\
                chrName=%r, name=%r, position=%r, weight=%r)' % \
                (repr(self.chrName), repr(self.name), repr(self.position),\
                repr(self.weight))

class DartTSSPair:
    """Object representing a dart-TSS pair
    
    A DartTSSPair is a Dart with an addition attribute of weight. This object
    is particularly useful when calculating the Gi local statistic.

    Parameters
    ----------
    chrName : str
              Name of the chromosome where the dart lives. (default = '')
    name : str
           Unique dart name. (default = '')
    position : int
               Unique dart position. (default = -1)
    weight : float
             Assigned weight. (default = -1)       

    Attributes
    ----------
    Same as parameters.

    Example
    --------
    """

    def __init__(self, chrName='', dartName='', dartPosition=-1,\
            TSSPosition=-1, weight=-1, geneName='', geneID=''):
        self.chrName = chrName
        self.dartName = dartName
        self.dartPosition = dartPosition
        self.TSSPosition = TSSPosition
        self.weight = weight
        self.geneName = geneName
        self.geneID = geneID

    def __str__(self):
        return "\t".join([self.chrName, self.dartName, str(self.dartPosition),\
                self.geneName, self.geneID, str(self.TSSPosition),\
                str(self.weight)])

    def __repr__(self):
        return 'GREATx.DartTSSPair(chrName=%r, dartName=%r, dartPosition=%r, \
                TSSPosition=%r, weight=%r, geneName=%r, geneID=%r)' %\
                (repr(self.chrName), repr(self.dartName),\
                repr(self.dartPosition), repr(self.TSSPosition),\
                repr(self.weight), repr(self.geneName), repr(self.geneID))

class TSS:
    """Object representing a TSS"""

    def __init__(self, position=-1, geneName='', geneID='', chrName=''):
        self.chrName = chrName
        self.geneName = geneName
        self.geneID = geneID
        self.position = int(position)

    def __str__(self):
        return 'chrName: %s\nGene Name: %s\nGene ID: %s\nPosition: %d' %\
                (self.chrName, self.geneName, self.geneID, self.position)

    def __repr__(self):
        return 'GREATx.TSS(position=%r, geneName=%r, geneID=%r, chrName=%r)'%\
                (repr(self.position), repr(self.geneName),\
                repr(self.geneID), repr(self.chrName))

class WeightedRegDom:
    """Object used to compute dart weights in a regulatory domain

    Weights are determined for each dart in the regulatory domain by overlapping
    normal distributions centered at each transcription start site. The
    regulatory domain is defined as all bases within the defined cut-off
    distance of a single transcription start site.

    Parameters
    ----------
    cutOff : int
             maximum distance an dart can be from
             a transcription start site to be
             assigned any weight from it
    mean : float
           mean value of the normal distribution
           for calculating weights
    sd : float
         standard deviation of the normal distribution
         for calculating weights

    Attributes
    ----------
    Same as parameters.

    Examples
    --------
    Find the best position for a dart in a regulatory domain
    
    >>> TSSs = [GREATx.TSS(chrName='chr123', geneName='tss', geneID='1234', position=i) for i in range(10)]
    >>> wgtRegDom = GREATx.WeightedRegDom(cutOff=2, mean=0.0, sd=3)
    >>> bestWeightedDart = wgtRegDom.bestWeightedDart(TSSs)
    
    Get a DartTSSPair
    
    >>> wgtRegDom = GREATx.WeightedRegDom(TSSs, cutOff=2, mean=0.0, sd=3)
    >>> Dart = Dart('myDart', 14)
    >>> TSS = TSS('chr123', 'myGene', 'myGeneID', 20)
    >>> DartTSSPair = wgtRegDom.makeDartTSSPair(Dart,TSS)
    """

    def __init__(self, cutOff=1000000, mean=0.0, sd=333333):
        self.cutOff = cutOff
        self.mean = mean
        self.sd = sd

    def __repr__(self):
        return 'WeightedRegDom(%r, %r, %r)' %\
                (repr(self.cutOff),repr(self.mean), repr(self.sd))

    # note that this algorithm is okay because transcriptionStartSites is sorted
    def bestWeightedDart(self, TSSs, chromosomes=HUMAN_CHROMOSOMES):
        """Returns the best WeightedDart for the given TSSs
        
        This method stupidly visits every base within the regulatory
        domain to determin which would be the best place to put a
        dart. A lot of work can be done to improve this algorithm.
        It is too slow to do anything reasonable with it.
        """

        bestWeightedDart = WeightedDart('', -1, weight=-1)
        for chrName in chromosomes:
            chrTSSs = filter(lambda x: x.chrName == chrName, TSSs)
            chrTSSs.sort(key=operator.attrgetter('position'))
            dart = Dart(chrName, -1)

            for tss in chrTSSs:
                for i in range(-self.cutOff, self.cutOff+1):

                    # dart position only increases
                    if (dart.position < tss.position + i):
                        dart.position = tss.position + i
                        wDart = self.getWeightedDart(chrTSSs, dart, wantFilter=False)

                        if (wDart.weight > bestWeightedDart.weight):
                            bestWeightedDart.chrName = wDart.chrName
                            bestWeightedDart.position = wDart.position
                            bestWeightedDart.weight = wDart.weight

        return bestWeightedDart

    def getWeightedDart(self, TSSs, dart, wantFilter=True):
        """Returns a WeightedDart for the given TSSs."""

        if wantFilter:
            TSSs = filter(lambda x: x.chrName == dart.chrName, TSSs)

        wDart = WeightedDart(dart.chrName, dart.position, weight=0)
        for tss in TSSs:
            if (fabs(wDart.position - tss.position) <= self.cutOff):
                wDart.weight += self.getDartTSSPairWgt(dart, tss)

        return wDart

    def getDartTSSPairWgt(self, dart, tss):
        """Return the weight for a particular dart relative to a TSS."""

        wgtDist = norm(self.mean, self.sd)
        maxWgtDist = wgtDist.pdf(self.mean)
        dartWgt = wgtDist.pdf(tss.position - dart.position)/maxWgtDist

        return dartWgt

    def makeDartTSSPair(self, dart, tss):
        if dart.chrName == tss.chrName:
            return DartTSSPair(chrName=dart.chrName,\
                    dartName=dart.name,\
                    dartPosition=dart.position,\
                    weight=self.getDartTSSPairWgt(dart, tss),\
                    geneName=tss.geneName,\
                    geneID=tss.geneID,\
                    TSSPosition=tss.position)
        else:
            sys.stderr('Cannot make a dart-TSS pair b/c \
                    dart.chrName != tss.chrName')
            return None

class TermDartTSSTriple:
    """Object used to parse data from the output file of AssociationMaker

    This object is primarily used to parse the output file of AssociationMaker
    into an array of objects, one for each line of the file. The name is
    awkward, but this object gets the job done.

    Parameters
    ----------
    termID : str
             unique term id
    chrName : str
              chromosome name
    dartName : str
               unique dart name
    dartPosition : int
                   unique dart position
    geneName : str
               unique gene name
    geneID : str
             unique gene id
    TSSPosition : int
                  unique TSS position
    weight : float
             weight for the dart-TSS pair
    percentCoverage : float
                      percent of the genome covered by this term id

    Attributes
    ----------
    Same as parameters.

    Example
    --------
    >>> lineObjects = []
    >>> for line in inFile:
    >>>     lineObjects.append(TermDartTSSTriple(line))
    >>> 
    >>> termIDs = list(set([lineObject.termID for lineObject in lineObjects]))
    >>> dartNames = list(set([lineObject.dartName for lineObject in lineObjects]))
    """

    def __init__(self, tripleLine):
        self._parseTermDartTSSTripleLine(tripleLine)

    def _parseTermDartTSSTripleLine(self, tripleLine):
        line = tripleLine.split()
        self.termID = line[0]
        self.chrName = line[1]
        self.dartName = line[2]
        self.dartPosition = int(line[3])
        self.geneName = line[4]
        self.geneID = line[5]
        self.TSSPosition = int(line[6])
        self.weight = float(line[7])
        self.percentCoverage = float(line[8])

    def __str__(self):
        return "\t".join([\
            self.termID,\
            self.chrName,\
            self.dartName,\
            str(self.dartPosition),\
            self.geneName,\
            self.geneID,\
            str(self.TSSPosition),\
            str(self.weight),\
            str(self.percentCoverage)])

    def __repr__(self):
        return "TermDataTSSTriple(%r)" % (repr("\t".join([\
                repr(self.termID),\
                repr(self.chrName),\
                repr(self.dartName),\
                repr(str(self.dartPosition)),\
                repr(self.geneName),\
                repr(self.geneID),\
                repr(str(self.TSSPosition)),\
                repr(self.weight),\
                repr(self.percentCoverage)])))


class Loci:
    """Object to represent a loci

    This object is a convenient way to represent a loci and its attributes.
    The primary role of this object is to parse the loci lines in a file
    like hg18.loci.

    Parameters
    ----------
    chrName : str
              chromosome name
    TSSPosition : int
                  unique TSS position
    geneName : str
               unique gene name
    geneID : str
             unique gene id
    strand : str
             +/- for upstream or down stream

    Attributes
    ----------
    Same as parameters.

    Example
    --------
    """

    def __init__(self, lociLine):
        self._parseLociLine(lociLine)

    def _parseLociLine(self, lociLine):
        """Expands fields of lociLine to a regulatory region line."""
        lineFields = lociLine.split()
        self.chrName = lineFields[1]
        self.TSSPosition = int(lineFields[2])
        self.geneName = lineFields[4]
        self.geneID = lineFields[0]
        self.strand = lineFields[3]

    def __str__(self):
        return "\t".join([\
                self.geneID,\
                self.chrName,\
                str(self.TSSPosition),\
                self.strand,\
                self.geneName])

    def __repr__(self):
        return "Loci(%r, %r)" % (repr("\t".join([\
                repr(self.geneID),\
                repr(self.chrName),\
                repr(self.TSSPosition),\
                repr(self.strand),\
                repr(self.geneName)])), repr(cutOff))

class LociRegulatoryRegion(Loci):
    """Object representing a loci expanded into a regulatory region

    This object uses its parent class to parse a loci line in a file
    like hg18.loci and can write a tab delimited line which includes
    the regulatory region for the loci specificed by cutOff in the
    paramaters during instantiation.

    Parameters
    ----------
    chrName : str
              chromosome name
    TSSPosition : int
                  unique TSS position
    geneName : str
               unique gene name
    geneID : str
             unique gene id
    strand : str
             +/- for upstream or down stream
    cutOff : int
            desired size of the regulatory region for each loci

    Attributes
    ----------
    Same as parameters.

    Example
    --------
    >>> loci = open(lociFn, 'r')
    >>> regDom = open(regDomFn, 'w')
    >>> 
    >>> for line in loci:
    >>>     regDom.write(str(LociRegulatoryRegion(line,cutOff=cutOff)) + '\n')
    """

    def __init__(self, *args, **kwargs):
        Loci.__init__(self, *args)
        cutOff = kwargs.pop('cutOff')
        self.regStart = max(0, int(self.TSSPosition) - cutOff)
        self.regEnd = int(self.TSSPosition) + cutOff

    def __str__(self):
        return "\t".join([\
                self.chrName,\
                str(self.regStart),\
                str(self.regEnd),\
                self.geneName,\
                self.geneID,\
                self.strand,\
                str(self.TSSPosition)])

    def __repr__(self):
        return "LociRegulatoryRegion(%r, cutOff=%r)" % (repr("\t".join([\
                repr(self.geneName),\
                repr(self.chrName),\
                repr(self.TSSPosition),\
                repr(self.strand),\
                repr(self.geneName),\
                repr(self.strand)])), repr(cutOff))

def createRegDomsFileFromTSSs(lociFn, regDomFn, cutOff):
    """Expands Loci file into a regulator regions file
    
    Parameters
    ----------
    
    lociFn : str
             name of file containing loci (e.g. hg18.loci)
    regDomFn : str
               name of file to contain the created regulatory regions
               (e.g. hg18.regDom.bed)
    cutOff : int
             cut-off for regulatory regions
    """
    loci = open(lociFn, 'r')
    regDom = open(regDomFn, 'w')
    
    for line in loci:
        regDom.write(str(LociRegulatoryRegion(line,cutOff=cutOff)) + '\n')

def overlapSelect(regDomFn, dartFn, mergedFn, options=''):
    """Runs overlapselect <options> <regDomFn> <dartFn> <mergedFn>.

    Parameters
    ----------
    selectFile : str
                 name of file containing loci (e.g. hg18.loci)
    inFile : str
             name of file to contain the created regulatory regions (e.g. hg18.regDom.bed)
    outFile : str
              cut-off for regulatory regions
    options : str
              space delemited list of options (e.g. "-selectFmt=bed -mergeOutput")

    Usage
    -----
    Note: This usage docstring is from the program overlapSelect
    
    overlapSelect [options] selectFile inFile outFile
    
    Select records based on overlaping chromosome ranges.  The ranges are
    specified in the selectFile, with each block specifying a range.
    Records are copied from the inFile to outFile based on the selection
    criteria.  Selection is based on blocks or exons rather than entire
    range.
    
    Options starting with -select* apply to selectFile and those starting
    with -in* apply to inFile.
    
    Options:
      -selectFmt=fmt - specify selectFile format:
              psl - PSL format (default for *.psl files).
              genePred - gepePred format (default for *.gp or
                         *.genePred files).
              bed - BED format (default for *.bed files).
                    If BED doesn't have blocks, the bed range is used. 
      -selectCoordCols=spec - selectFile is tab-separate with coordinates
           as described by spec, which is one of:
                o chromCol - chrom in this column followed by start and end.
                o chromCol,startCol,endCol - chrom, start, and end in specified
                  columns.
                o chromCol,startCol,endCol,strandCol - chrom, start, end, and
                  strand in specified columns.
              NOTE: column numbers are zero-based
      -selectCds - Use only CDS in the selectFile
      -selectRange - Use entire range instead of blocks from records in
              the selectFile.
      -inFmt=fmt - specify inFile format, same values as -selectFmt.
      -inCoordCols=spec - inFile is tab-separate with coordinates specified by
          spec, in format described above.
      -inCds - Use only CDS in the inFile
      -inRange - Use entire range instead of blocks of records in the inFile.
      -nonOverlapping - select non-overlaping instead of overlaping records
      -strand - must be on the same strand to be considered overlaping
      -oppositeStrand - must be on the opposite strand to be considered overlaping
      -excludeSelf - don't compare records with the same coordinates and name.
          Warning: using only one of -inCds or -selectCds will result in different
          coordinates for the same record.
      -idMatch - only select overlapping records if they have the same id
      -aggregate - instead of computing overlap bases on individual select entries, 
          compute it based on the total number of inFile bases overlap by selectFile
          records. -overlapSimilarity and -mergeOutput will not work with
          this option.
      -overlapThreshold=0.0 - minimun fraction of an inFile record that
          must be overlapped by a single select record to be considered
          overlapping.  Note that this is only coverage by a single select
          record, not total coverage.
      -overlapThresholdCeil=1.1 - select only inFile records with less than
          this amount of overlap with a single record.
      -overlapSimilarity=0.0 - minimun fraction of inFile and select records that
          Note that this is only coverage by a single select record and this
          is; bidirectional inFile and selectFile must overlap by this
          amount.  A value of 1.0 will select identical records (or CDS if
          both CDS options are specified.  Not currently supported with
          -aggregate.
      -overlapSimilarityCeil=1.1 - select only inFile records with less than
          this amount of similarity with a single record.
      -overlapBases=-1 - minimum number of bases of overlap, < 0 disables.
      -statsOutput - output overlap statistics instead of selected records. 
          If no overlap criteria is specified, all overlapping entries are
          reported, Otherwise only the pairs passing the citeria are
          reported. This results in a tab-seperated file with the columns:
             inId selectId inOverlap selectOverlap overBases
          Where inOverlap is the fraction of the inFile record overlapped by
          the selectFile record and selectOverlap is the fraction of the
          select record overlap by inFile records.  With -aggregate, output
          is:
             inId inOverlap inOverBases inBases
      -statsOutputAll - like -statsOutput, however output all inFile records,
          even ones that are not overlapped.
      -mergeOutput - output file with be a merge of the input file with the
          selectFile records that selected it.  The format is
             inRec<tab>selectRec.
          if multiple select records hit, inRec is repeated. This will increase
          the memory required. Not supported with -nonOverlapping or -aggregate.
      -idOutput - output a table seprate file of pairs of
             inId selectId
          with -aggregate, omly a single column of inId is written
      -dropped=file  - output rows that were dropped to this file.
      -verbose=n - verbose > 1 prints some details,
    """
    
    program = "/afs/ir/class/cs173/bin/i386_linux26/overlapSelect"
    os.system(" ".join([program, options, regDomFn, dartFn, mergedFn]))

def assignWeights(cutOff, mean, sd, mergedFn, dartsToWeightsFn):
    """Writes to dartsToWeightsFn each dart with the geneName, geneID, and weight

     merged file must follow this format:
    |<-----------Dart Information---------->|<---------RegDom Info------->| Gene Name | GeneID | strand | TSS Location |
    | chrName | chrStart | chrEnd | dartName | chrName | chrStart | chrEnd | Gene Name | GeneID | strand | TSS Location |

    dartsToWeights file will be in this format:
    | dartName | Gene Name | Gene ID | Weight [0-1] |

    """

    merged = open(mergedFn, 'r')
    dartsToWeightsFile = open(dartsToWeightsFn, 'w')
    wgtRegDom = WeightedRegDom(cutOff, mean, sd)

    dart = Dart()
    tss = TSS()
    for line in merged:
        line = line.split()

        dart.chrName = line[0]
        dart.name = line[3] 
        dart.position = (int(line[1]) + int(line[2]))/2

        tss.chrName = line[4]
        tss.position = int(line[10])
        tss.geneName = line[7]
        tss.geneID = line[8]

        dartTSSPair = wgtRegDom.makeDartTSSPair(dart, tss)
        dartsToWeightsFile.write(str(dartTSSPair) + "\n")


class RegDom:
    """Instantiates a regulatory domain"""
    def __init__(self, start, end, id):
        self.id = id
        self.start = start
        self.end = end


class AssociationMaker:
    # Output format:
    """| term id# |  chrom name | arrow | arrow position (relative to chrom) | gene name | gene id | dart weight |  term coverage % (between 0 and 1) |"""
    dartTSSPairs = []
    genetoterms = collections.defaultdict(lambda :[])
    termtocoverage = collections.defaultdict(lambda : 0.0)

    def __init__(self, dartsToWeightsFn, geneOntologyFn, regDomFn):
        self.readDartWeightsFile(dartsToWeightsFn)
        self.buildGeneTermMap(geneOntologyFn)
        self.buildTermWeightsMap(regDomFn)

    def readDartWeightsFile(self, fstr):
        with open(fstr) as f:
            for line in f:
                line = line.split("\t")
                dartTSSPair = DartTSSPair(\
                        chrName=line[0],\
                        dartName=line[1],\
                        dartPosition=line[2],\
                        geneName=line[3],\
                        geneID=line[4],\
                        TSSPosition=line[5],\
                        weight=float(line[6]))

                self.dartTSSPairs.append(dartTSSPair)

    def buildGeneTermMap(self, geneOntologyFn):
        f = open(geneOntologyFn)
        for line in f:
            line = line.split("\t")
            term_id = int(re.search("\d+", line[0]).group(0))
            geneID = re.search("\d+", line[1]).group(0)
            self.genetoterms[geneID].append(term_id)

    def buildTermWeightsMap(self, regDomFn):
        genes = [dartTSSPair.geneID for dartTSSPair in self.dartTSSPairs]
        regdoms = []
        with open(regDomFn) as f:
            for line in f:
                line = line.split()
                if line[4] in genes:
                    regdoms.append(RegDom(int(line[1]),int(line[2]),line[4]))
        genome_size = sum(HUMAN_CHROMOSOME_SIZES)
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

    def getTerms(self, geneID):
        #return [str(x[0]) for x in self.termtogenes.items() if gene in x[1]]
        return map(lambda x: str(x), self.genetoterms[geneID])

    def buildLine(self, term_name, dartTSSPair, term_weight):
        return term_name + "\t" + str(dartTSSPair) + "\t" + str(term_weight) + "\n"

    def writeOutput(self, output_file):
        f = open(output_file, "w")
        for dartTSSPair in self.dartTSSPairs:
            terms = self.getTerms(dartTSSPair.geneID)
            if terms == []: # in case we get a gene that for some reason has no terms associated
                f.write(self.buildLine("UNKNOWN", dartTSSPair, "0.0"))
            else:
                for term in terms:
                    f.write(self.buildLine(term, dartTSSPair, self.termtocoverage[term]))
        f.close()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="%prog <lociFn> <ontoToGeneFn> <dartFn> <SRFtoTermsFn> <outFn> \
            <cutOff> <mean> <sd> <which beta>",
                          description=(""))

    #parser.add_option("-q", "--quiet",
    #                  action="store_false", dest="verbose", default=True,
    #                  help="don't print status messages to stdout")
    
    (options, args) = parser.parse_args()
    if (len(args) != 6):
        parser.print_usage()
        sys.exit(1)

    lociFn = args[0]
    ontoToGeneFn = args[1]
    dartFn = args[2]
    SRFtoTermsFn = args[3]
    outFn = args[4]
    cutOff = int(args[5])
    mean = float(args[6])
    sd = float(args[7])
    whichBeta = int(args[8])

    # get data/SRFtoTerms.data
    createRegDomsFileFromTSSs(lociFn, "/tmp/hg18.regDom.bed", 1000000)
    overlapSelect("/tmp/hg18.regDom.bed", dartFn, "/tmp/regDom.SRF.merge", options="-mergeOutput")
    assignWeights(cutOff, mean, sd, "/tmp/regDom.SRF.merge", "/tmp/SRF.wgt")
    maker = AssociationMaker("/tmp/SRF.wgt", ontoToGeneFn, "/tmp/hg18.regDom.bed")
    maker.writeOutput(SRFtoTermsFn)

    # remove /tmp files
    os.system("rm /tmp/hg18.regDom.bed /tmp/regDom.SRF.merge /tmp/SRF.wgt")
    
    inFile = open(SRFtoTermsFn)
    outFile = open(outFn, 'w')
    
    weightSum = 0 #Will store our alpha
    X = 0 #Will store our input probability 
    
    lineObjects = []
    for line in inFile:
        lineObjects.append(TermDartTSSTriple(line))
    
    termIDs = list(set([lineObject.termID for lineObject in lineObjects]))
    dartNames = list(set([lineObject.dartName for lineObject in lineObjects]))
    
    termPValues = {}
    
    for termID in termIDs:
        if termID != 'UNKNOWN':
            termIDObjects = filter(lambda x: x.termID == termID, lineObjects)
            weights = [termIDObject.weight for termIDObject in termIDObjects]
            genes = list(set([termIDObject.geneName for termIDObject in termIDObjects]))
    
        
            alpha = sum(weights)
            x = termIDObjects[0].percentCoverage
            
            ## new possibility for x?
            #wgtDist = norm(self.mean, self.sd)
            #x = len(genes)*len(wgtDist.cdf(cutOff) - wgtDist.cdf(-cutOff))\
            #        /sum(HUMAN_CHROMOSOME_SIZES)
            
            if whichBeta == 1:
                beta = len(termIDObjects) - alpha
            elif whichBeta == 2:
                beta = len(termIDObjects) * max(weights) - alpha
            elif whichBeta == 3:
                wgtRegDom = WeightedRegDom(cutOff=1000000, mean=0, sd=333333)
                termTSSs = [TSS(position=termIDObject.TSSPosition,\
                                       geneName=termIDObject.geneName,\
                                       geneID=termIDObject.geneID,\
                                       chrName=termIDObject.chrName)\
                            for termIDObject in termIDObjects]
                beta = len(termIDObjects) * (wgtRegDom.bestWeightedDart(termTSSs, chromosomes=HUMAN_CHROMOSOMES)).weight
            elif whichBeta == 4:
                for dartName in dartNames:
                    dartNameObjects = filter(lambda x: x.dartName == dartName, termIDObjets)
                    weight = sum([dartNameObject.weight for dartNameObject in dartNameObjects]) 
                    weightedDart = GREATx.WeightedDart(chrName=dartNameObjects[0].chrName,\
                                                       name=dartName,
                                                       position=dartNameObjects[0].dartPosition,
                                                       weight=weight)
                    weightedDarts += weightedDart

                beta = len(weightedDarts) * max([weightedDart.weight for weightedDart in weightedDarts]) - alpha

    
            outFile.write(termID + "\t" + str(scipy.stats.beta.cdf(x, alpha, beta)) + "\n")

