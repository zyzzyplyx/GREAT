import sys
import os
import WeightedRegDom as wrd #Our own module that computes the maximum weight of a given set of transcription start sites

inFileName = sys.argv[1]
outFileName = sys.argv[2]
inFile = open(inFileName)
outFile = open(outFileName,'w')
#outFile after completion of this script will contain [X \alpha \beta] on one line



weightSum = 0 #Will store our alpha
X = 0 #Will store our input probability 
weightIndex = 3 #the column with the weight
percentIndex = 4 #the column with the percent
tssIndex = 5#the column with the TSS

for line in inFile:
	line = line.split() #Is it tab delimited??
	weightSum = float(line[weightIndex])

X = float(line[percentIndex])

###inFile will have the transcription start sites for the genes in the last column
#I need to do a sort then uniq on the transcription start site, this final file will be
#called tss


os.system("awk '{print $6}'" + inFile "> temp.txt") #cuts the transcription start sites
os.system("sort -g -u temp.txt > tss.txt") #Sorts them then uniques


#Open the transcription start site file
tssFile = open("tss.txt",'w')
 
#TO DO: Update this sequence to be the correct
transcriptionStartSites = wrd.readPositions(tssFile)
wgtRegDom = wrd.WeightedRegDom(transcriptionStartSites, cutoff=CUTOFF, mean=MEAN, sd=SD)
r = wgtRegDom.maxArrowWgt()
N = wrd.repr(r.maxWgt) 

outFile.write("%f \t %f \t %f\n",X,weightSum, N)
inFile.close()
tssFile.close()
outFile.close()

#Delete intermediary files
os.system("rm -f temp.txt")
os.system("rm -f tss.txt")



###
## python commands with new GREATx module

#import GREATx
#GREATx.createRegDomsFileFromTSSs("/afs/ir/class/cs173/finalProjects/GREATRegDoms/ontologies/hg18/hg18.loci", "/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/hg18.regDom.bed", 1000000)
#GREATx.overlapSelect("/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/hg18.regDom.bed", "/afs/ir/class/cs173/finalProjects/GREATRegDoms/SRF.hg18.bed", "/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/regDom.SRF.merge", options="-mergeOutput")
#GREATx.assignWeights(1000000, 0, 333333, "/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/regDom.SRF.merge", "/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/SRF.wgt")
#maker = GREATx.AssociationMaker("/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/SRF.wgt", "/afs/ir/class/cs173/finalProjects/GREATRegDoms/ontologies/hg18/GOBiologicalProcess/ontoToGene.canon", "/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/hg18.regDom.bed", "/afs/ir/class/cs173/finalProjects/GREATRegDoms/ontologies/hg18/chrom.sizes")
#maker.writeOutput("/afs/ir.stanford.edu/users/c/c/cceleri/GREAT/data/SRFtoTerms.data")


