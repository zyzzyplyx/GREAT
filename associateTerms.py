#!/usr/bin/env python

import collections
import re
import sys
#Usage: ./associateTerms.py dart_weights_file output_file [gene to term files]

class Dart:

    def __init__(self, name, gene_name, id, weight):
        self.name = name
        self.gene_name = gene_name
        self.gene_id = id
        self.weight = weight

    def __str__(self):
        return 'Dart Name: ' + self.name + '\nGene Name: ' + self.gene_name + '\nID: ' + str(self.id) + '\nWeight: ' + str(self.weight)

class AssociationMaker:

    darts = []
    genetoterms = collections.defaultdict(lambda :[])

    def __init__(self, filestr, gene_term_files):
        self.readDartWeightsFile(filestr)
        self.buildGeneTermMap(gene_term_files)

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

    def getTerms(self, gene_id):
        #return [str(x[0]) for x in self.termtogenes.items() if gene in x[1]]
        return map(lambda x: str(x), self.genetoterms[gene_id])

    def buildLine(self, term_name, dart):
        return "\t".join([term_name, dart.gene_name, dart.name, str(dart.weight)]) + "\n"

    def writeOutput(self, output_file):
        f = open(output_file, "w")
        for dart in self.darts:
            terms = self.getTerms(dart.gene_id)
            if terms == []: # in case we get a gene that for some reason has no terms associated
                f.write(self.buildLine("UNKNOWN", dart))
            else:
                for term in terms:
                    f.write(self.buildLine(term, dart))
        f.close()

def main():
   maker = AssociationMaker(sys.argv[1],sys.argv[3:])
   maker.writeOutput(sys.argv[2])

if __name__ == '__main__':
    main()
