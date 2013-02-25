This directory holds the GREAT-related code and makefile required to create working
copies of the novel sections of the core calculation engine functionality of
GREAT (McLean CY et al., Nature Biotechnology, 2010).
These functions include the following steps:
 1. Computational calculation of the regulatory domain of all genes
 2. Calculating a genomic region-based binomial p-value given a set of
regulatory domains associated with the term.

The GREAT codebase is dependent on the UCSC Kent source tree.  Consequently, the
UCSC Kent source tree must be downloaded and installed before attempting to build the GREAT
core calculation engine tools.  All steps required to get the GREAT core calculation engine
tools up and running are documented below.  These instructions are intended for Unix/Linux
systems only; building the GREAT core engine tools on other systems is beyond the scope of
these instructions.


### Building the provided GREAT core calculation engine tools ###
1. Download and install the UCSC Kent source tree
     a. A current copy of the Kent source tree is freely available for academic, nonprofit,
           and personal use at http://hgdownload.cse.ucsc.edu/admin/jksrc.zip, and can also
           be obtained via CVS (see http://genome.ucsc.edu/admin/cvs.html).
     b. Build and installation instructions for the Kent source tree are available at
           http://genome.ucsc.edu/admin/jk-install.html.  Only the first six steps are required;
           in fact only the library jkweb.a is necessary.

2. Update the GREAT core calculation engine tools makefile to access the required Kent source library functions
     a. From step 1, you have installed the Kent source tree in some base directory (hereafter
           referenced as $BASE_DIR).  Within the $BASE_DIR/kent/src/lib/$MACHTYPE/ directory
           there should be a library file named jkweb.a.  Verify that this file exists.  If it
           does not exist, step 1 was not performed properly and should be redone.
     b. Open the GREAT core calculation engine tools makefile using a text editor.
		   The first line reads "KENT_DIR = path/to/your/kent/src".
           Update this path definition with the actual location of your Kent source directory.  This will be
           $BASE_DIR/kent/src for whatever value of $BASE_DIR is appropriate.

3. Build the GREAT core calculation engine tools
     a. If the KENT_DIR assignment within the makefile is set up correctly, simply typing 'make' in the
           GREAT directory will build the GREAT core calculation engine tools.
			* An executable named createRegulatoryDomains can be used to calculate regulatory domains.
			* An executable named calculateBinomialP can be used to calculate genomic region-based binomial p-values
           Calling either program with no arguments prints its usage message.



### Running the createRegulatoryDomains tool ###
The createRegulatoryDomains tool is used to generate computationally defined regulatory domains for all genes in a gene set.
The tool requires four arguments:

1. TSS.in
	This is a file holding a list of all genes to which you want to assign regulatory domains.  Each line of the file
	should correspond to a single gene to which you assign a regulatory domain.  Each line should have four fields
	tab-delimited:
              chromosome      transcription start site      strand      geneName

2. chrom.sizes
	This is a tab-delimited file holding the number of basepairs in each chromosome, in the following two-field format:
              chromosome      chromosome size

3. oneClosest|twoClosest|basalPlusExtension
	This argument corresponds to the type of association rule desired (see the
	"Association rules from genomic regions to genes" section of the Online Methods of our paper at
	 http://dx.doi.org/10.1038/nbt.1630 for a full description of each method).

4. regDoms.out
	This is the output file listing all of the computationally defined regulatory domains of each gene.  This file is
	in the following tab-delimited format:

              chromosome      chromStart      chromEnd      geneName     transcription start site      strand

	The chromosome, geneName, transcription start site, and strand are all identical to those in the input TSS.in file.
	The span of [chromStart, chromEnd) is the computationally-defined regulatory domain of the gene.  Note that this
	file is a valid BED file (http://genome.ucsc.edu/FAQ/FAQformat.html#format1).


Options:
	The -maxExtension, -basalUpstream, and -basalDownstream options allow users to vary the amount of genome
	associated with each gene.  Note that the basalUpstream and basalDownstream options are only relevant to the
	basalPlusExtension association rule.


### Running the calculateBinomialP tool ###
The calculateBinomialP tool is used to calculate the genomic region-based binomial p-value of enrichment for a particular
ontology term, based on the fraction of the genome associated with the term, the number of genomic regions in the input
set, and the number of genomic regions that are associated with genes annotated with the term.  The tool requires four
arguments:

1. regdoms.in
	This is a file of fully-specified regulatory domains of all genes in the genome that are annotated with the term
	of interest.  The format of the file follows the output of the createRegulatoryDomains tool:

              chromosome      chromStart      chromEnd      geneName     transcription start site      strand

	Note that these regulatory domains may overlap each other.

2. antigap.bed
	This is a BED file (http://genome.ucsc.edu/FAQ/FAQformat.html#format1) specifying all regions of the genome in which
	input genomic regions may land (e.g. all non–assembly gap base pairs in the genome).

	Important note:  antigap.bed is required to consist entirely of non-overlapping regions.

3. numTotalRegions
	The total number of genomic regions in the input set.

4. numRegionsHit
	The number of input genomic regions that are annotated with the ontology term of interest (due to their midpoints
	overlapping the regulatory domain of one or more genes annotated with the term).


The binomial p-value of enrichment for the term, given the four inputs, is printed to standard output.


### Additional help ###
Please direct any questions about compilation or usage to great@cs.stanford.edu.
