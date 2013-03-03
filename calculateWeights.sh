#!/bin/bash

#From a list of transcription start sites, creates regions of +/- 1MB
# Output is in BED format.
./createRegDoms.py GREATRegDoms/ontologies/hg18/hg18.loci hg18_regions.bed

# Using the regions output in the last step, identifies all overlaps between an
# arrow and a regdom.  One arrow may overlap with more than one regdom.
# Outfile Format:
#|<--------------Arrow Information--------------->|<-------RegDom Info---------->| Gene Name |  GeneID | strand | TSS Location |
#| chr5    60663738        60663788        SRF.1  | chr5   59031713    61031713  |  DEPDC1B  |  19474  |    -   |  60031713    |
 /afs/ir/class/cs173/bin/i386_linux26/overlapSelect -mergeOutput hg18_regions.bed GREATRegDoms/SRF.hg18.bed try.SRF.overlap.regDom.txt

# For each overlap, calculates the "weight" of the arrow, based on how close
# a given arrow is to the TSS.
# Outfile Format:
# | “Arrow” name (e.g. SRF.1) | Gene Name | Gene ID (numerical ID) | Weight (in the range of 0 to 1) |
# |           SRF.1           |  DEPDC1B  |        19474           |          0.165703737676         |
./findDist.py try.SRF.overlap.regDom.txt hg18_weights.bed

