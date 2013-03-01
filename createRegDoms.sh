#!/bin/bash

#From a list of transcription start sites, creates regions of +/- 1MB
cat GREATRegDoms/ontologies/hg18/hg18.loci | awk 'BEGIN{} {print $2,$3-1000000,$3+1000000, $4,$5}' > hg18_regions.bed



