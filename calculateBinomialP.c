#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "linefile.h"
#include "basicBed.h"
#include "genomeRangeTree.h"
#include "options.h"
#include "regdom.h"

#define MAXIT 10000
#define EPS 3.0e-7 
#define FPMIN 1.0e-30 

struct optionSpec options[] = {
	{NULL, 0}
};

void usage()
/* Explain usage and exit. */
/* example: createRegulatoryDomains TSS.in chrom.sizes basalPlusExtension regDoms.out */
{
  errAbort(
  "\n"
  "Calculates the binomial p-value of enrichment for a term given the list of regulatory domains of genes associated\n"
  "with the term, a list of all valid regions in the genome to include in the weight calculation of the term,\n"
  "the number of genomic regions in the entire input set, and the number of genomic regions that hit one of the input\n"
  "regulatory domains.  P-value is printed to standard output.\n\n"
  "Usage:\n\n"
  "calculateBinomialP regdoms.in antigap.bed numTotalRegions numRegionsHit\n"
  );
}

struct genomeRangeTree* getRangeTreeOfRegdoms(struct regdom* regdoms)
{
	struct genomeRangeTree *ranges = genomeRangeTreeNew();
	struct regdom* currRD;
	for (currRD = regdoms; currRD != NULL; currRD = currRD->next) {
		genomeRangeTreeAdd(ranges, currRD->chrom, currRD->chromStart, currRD->chromEnd);
	}
	return ranges;
}

long getTotalNonGapBases(struct bed* antigaps)
{
	long retval = 0;
	struct bed* curr;
	for (curr = antigaps; curr != NULL; curr = curr->next) {
		retval += (curr->chromEnd - curr->chromStart);
	}
	return retval;
}

long getAnnotatedNonGapBases(struct genomeRangeTree* ranges, struct bed* antigaps)
{
	long retval = 0;
	struct bed* currAntigap;
	for (currAntigap = antigaps; currAntigap != NULL; currAntigap = currAntigap->next) {
		retval += genomeRangeTreeOverlapSize(ranges, currAntigap->chrom, currAntigap->chromStart, currAntigap->chromEnd);
	}
	return retval;
}


double betacf(double a, double b, double x)
// Used by betai: Evaluates continued fraction for incomplete beta function 
// by modified Lentz s method (§5.2). 
{
	int m,m2;
	double aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	// These q s will be used in factors that occur in the coefficients (6.4.6). 
	qap=a+1.0;
	qam=a-1.0;
	c=1.0; // First step of Lentz s method. 
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN)
		d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1; m<=MAXIT; m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d; // One step (the even one) of the recurrence. 
		if (fabs(d) < FPMIN)
			d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN)
			c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d; // Next step of the recurrence (the odd one). 
		if (fabs(d) < FPMIN)
			d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN)
			c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS)
			break;
	}
	if (m > MAXIT)
		errAbort("a or b too big, or MAXIT too small in betacf");
	return h;
}


double betai(double a, double b, double x)
// Returns the incomplete beta function Ix(a, b). 
{
	double bt;

	if (x < 0.0 || x > 1.0)
		errAbort("Bad x in routine betai");
	if (x == 0.0 || x == 1.0)
		bt=0.0;
	else // Factors in front of the continued fraction. 
		bt=exp(lgamma(a+b)-lgamma(a)-lgamma(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0)) // Use continued fraction directly. 
		return bt*betacf(a,b,x)/a;
	else // Use continued fraction after making the symmetry transformation. 
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}


double getBinomPval(int n, int k, double p)
{
	if (k == 0) return 1;
	else return betai(k, n-k+1, p);
}

void calculateBinomialP(char* regdomFn, char* antigapFn, int totalRegions, int hitRegions)
/* Calculate binomial p-value of enrichment based on regulatory domains and regions hit */
{
	struct regdom* regdoms = readInitializedRegdomFile(regdomFn);

	// This will hold the union of all regulatory domains for quick search
	struct genomeRangeTree *ranges = getRangeTreeOfRegdoms(regdoms);

	// NOTE: Each of these regions must be non-overlapping.
	struct bed* antigaps = bedLoadAll(antigapFn);
	long totalNonGapBases = getTotalNonGapBases(antigaps);
	long annotatedNonGapBases = getAnnotatedNonGapBases(ranges, antigaps);

	double annotationWeight = (double)annotatedNonGapBases/(double)totalNonGapBases;

	double binomP = getBinomPval(totalRegions, hitRegions, annotationWeight);

	printf("%e\n", binomP);

	regdomFreeList(&regdoms);
	bedFreeList(&antigaps);
	genomeRangeTreeFree(&ranges);
}


int main(int argc, char *argv[])
{
	/* Get all inputs */
	optionInit(&argc, argv, options);
	if (argc != 5) usage();

	char *regdomFn, *antigapFn;
	int totalRegions, hitRegions;
	regdomFn = argv[1];
	antigapFn = argv[2];
	totalRegions = intExp(argv[3]);
	hitRegions = intExp(argv[4]);

	calculateBinomialP(regdomFn, antigapFn, totalRegions, hitRegions);

	optionFree();

	return 0;
}

