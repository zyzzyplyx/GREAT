#!/bin/bash

#From a list of transcription start sites, creates regions of +/- 1MB
./createRegDoms.py GREATRegDoms/ontologies/hg18/hg18.loci hg18_regions.bed
 /afs/ir/class/cs173/bin/i386_linux26/overlapSelect -mergeOutput hg18_regions.bed GREATRegDoms/SRF.hg18.bed try.SRF.overlap.regDom.txt
./findDist.py try.SRF.overlap.regDom.txt hg18_weights.bed

