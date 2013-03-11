import sys
import os
import GREATx 

# get data/SRFtoTerms.data
GREATx.createRegDomsFileFromTSSs("GREATRegDoms/ontologies/hg18/hg18.loci", "data/hg18.regDom.bed", 1000000)
GREATx.overlapSelect("data/hg18.regDom.bed", "GREATRegDoms/SRF.hg18.bed", "data/regDom.SRF.merge", options="-mergeOutput")
GREATx.assignWeights(1000000, 0, 333333, "data/regDom.SRF.merge", "data/SRF.wgt")
maker = GREATx.AssociationMaker("data/SRF.wgt", "GREATRegDoms/ontologies/hg18/GOBiologicalProcess/ontoToGene.canon", "data/hg18.regDom.bed", "GREATRegDoms/ontologies/hg18/chrom.sizes")
maker.writeOutput("data/SRFtoTerms.data")

inFile = open("data/SRFtoTerms.data")
outFile = open("data/Term.p_values", 'w')

weightSum = 0 #Will store our alpha
X = 0 #Will store our input probability 

weightIndex = 7 #the column with the weight
percentIndex = 8 #the column with the percent
tssIndex = 6 #the column with the TSS

for line in inFile:
	line = line.split() #Is it tab delimited??
	weightSum += float(line[weightIndex]) 

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


