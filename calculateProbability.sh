#!/bin/bash

#Usage: ./calculateProbability.sh inFile outFile
#Infile is assumed to have weights, %covered by genome, transcription start sites
#Outfile will ouput the p-value

inFile=$1
outFile=$2



python calculateStats.py inFile temp

export X=$(awk '{print $1}' temp)
export alpha=$(awk '{print $1}' temp)
export beta=$(awk '{print $3}' temp)

rm-f temp

./betaCDF $X $alpha $beta > outFile
