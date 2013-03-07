import sys
import os
import maxWeight

inFileName = sys.argv[1]
outFileName = sys.argv[2]
inFile = open(inFileName)
outFile = open(outFileName,'w')
#outFile will ouput [X \alpha \beta] on one line

tssFile = open("tempFile",'w')

weightSum = 0
X = 0
weightIndex = 0 #the column with the weight
percentIndex = 4 #the column with the percent

for line in inFile:
	line = line.split() #Is it tab delimited??
	weightSum = float(weightSum) + line[weightIndex]

X = line[percentIndex]

#########################################
#####INSERT CODE FROM MAXWEIGHT HERE#####
#N should store the maximum weight for the regions

# N = maxWeight.getWeights() ??

#######################################
#######################################



outFile.write("%f \t %f \t %f\n",X,weightSum, N)
inFile.close()
tssFile.close()
