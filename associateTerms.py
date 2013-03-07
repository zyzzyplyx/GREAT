#!/usr/bin/env python

import collections
import re
import sys
#Usage: ./associateTerms.py dart_weights_file output_file [gene to term files]

# Output format:
# | term id# |  gene name | arrow |      weight    |  term coverage % (between 0 and 1) |
# |   7154   |   DEPDC1B  | SRF.1 | 0.165703737676 |           0.00257426987979         |

class Dart:

    def __init__(self, name, gene_name, id, weight):
        self.name = name
        self.gene_name = gene_name
        self.gene_id = id
        self.weight = weight

    def __str__(self):
        return 'Dart Name: ' + self.name + '\nGene Name: ' + self.gene_name + '\nID: ' + str(self.id) + '\nWeight: ' + str(self.weight)

class RegDom:

    def __init__(self, start, end, id):
        self.id = id
        self.start = start
        self.end = end

class AssociationMaker:

    darts = []
    genetoterms = collections.defaultdict(lambda :[])
    termtocoverage = collections.defaultdict(lambda : 0.0)

    def __init__(self, filestr, gene_term_files):
        self.readDartWeightsFile(filestr)
        self.buildGeneTermMap(gene_term_files)
        self.buildTermWeightsMap()

    def readDartWeightsFile(self, fstr):
        with open(fstr) as f:
            for line in f:
                line = line.split("\t")
                dart = Dart(line[0], line[1], int(line[2]), float(line[3]))
                self.darts.append(dart)

    def buildGeneTermMap(self, gene_term_files):
        for mfile in gene_term_files:
            with open(mfile) as f:
                for line in f:
                    line = line.split("\t")
                    term_id = int(re.search("\d+", line[0]).group(0))
                    gene_id = int(re.search("\d+", line[1]).group(0))
                    self.genetoterms[gene_id].append(term_id)

    def buildTermWeightsMap(self):
        genes = [dart.gene_id for dart in self.darts]
        regdoms = []
        with open("hg18_regions.bed") as f:
            for line in f:
                line = line.split()
                if int(line[4]) in genes:
                    regdoms.append(RegDom(int(line[1]),int(line[2]),int(line[4])))
        genome_size = 0
        with open("GREATRegDoms/ontologies/hg18/chrom.sizes") as f:
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

def main():
   maker = AssociationMaker(sys.argv[1],sys.argv[3:])
   maker.writeOutput(sys.argv[2])

if __name__ == '__main__':
    main()
