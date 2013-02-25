#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "regdom.h"

struct optionSpec options[] = {
	{"maxExtension", OPTION_INT},
	{"basalUpstream", OPTION_INT},
	{"basalDownstream", OPTION_INT},
	{NULL, 0}
};

int maxExtension = 1000000;   // Maximum distance to extend regulatory domain in absence of other genes
int basalUpstream = 5000;     // Basal regulatory domain extension upstream of TSS
int basalDownstream = 1000;   // Basal regulatory domain extension downstream of TSS

// The three association rule types
const char* oneClosest = "oneClosest";
const char* twoClosest = "twoClosest";
const char* basalPlusExtension = "basalPlusExtension";

struct hash* chromSizesHash;

void usage()
/* Explain usage and exit. */
/* example: createRegulatoryDomains TSS.in chrom.sizes basalPlusExtension regDoms.out */
{
  errAbort(
  "\n"
  "Creates regulatory domains for a set of genes based on the genomic location of the TSS and the association rule used.\n\n"
  "Usage:\n\n"
  "createRegulatoryDomains TSS.in chrom.sizes [%s|%s|%s] regDoms.out [options]\n\n"
  "Options:\n\n"
  "-maxExtension=N	 Distance to extend a gene's regulatory region from the TSS in absence of any other nearby genes (default=%d)\n"
  "-basalUpstream=N	  Basal regulatory region extension distance upstream (strand-dependent!) of TSS (default=%d)\n"
  "-basalDownstream=N	Basal regulatory region extension distance downstream (strand-dependent!) of TSS (default=%d)\n",
  oneClosest, twoClosest, basalPlusExtension, maxExtension, basalUpstream, basalDownstream
  );
}

void validateInput(char* association)
{
	if (!sameString(association, oneClosest) && !sameString(association, twoClosest) && !sameString(association, basalPlusExtension))
		errAbort("Association rule must be one of %s, %s, %s", oneClosest, twoClosest, basalPlusExtension);
	if (maxExtension < 0)
		errAbort("Maximum extension must be a non-negative integer: %d", maxExtension);
	if (basalUpstream < 0)
		errAbort("Basal upstream must be a non-negative integer: %d", basalUpstream);
	if (basalDownstream < 0)
		errAbort("Basal downstream must be a non-negative integer: %d", basalDownstream);

	if ((sameString(association, oneClosest) || sameString(association, twoClosest)) &&
	   (optionExists("basalUpstream") || optionExists("basalDownstream")))
		errAbort("Basal up/downstream options only relevant to basalPlusExtension association rule");
}

void readChromSizes(char* fn)
/* Reads a file of "chrom \t chromSize" elements into global hash keyed by chrom */
{
	int fieldCount = 2;
	char* row[fieldCount];
	chromSizesHash = newHash(5);  // Expect 2^5 = 32 or fewer chromosomes
	int wordCount;
	int chromSize = -1;
	struct lineFile* lf = lineFileOpen(fn, TRUE);
	while ((wordCount = lineFileChopTab(lf, row)) != 0) {
		if (wordCount != fieldCount)
			errAbort("Expecting exactly %d words line %d of %s", fieldCount, lf->lineIx, lf->fileName);
		chromSize = lineFileNeedNum(lf, row, 1);
		hashAddInt(chromSizesHash, row[0], chromSize);
	}
	lineFileClose(&lf);
}

void createBasalPlusExtensionRegDoms(struct regdom* regdoms, int maximumExtension, int basalUp, int basalDown)
{
	struct regdom *prev = NULL, *curr = NULL, *next = NULL;
	int currChromSize = -1;

	// First map all regdoms to their basal regions
	for (curr = regdoms; curr != NULL; curr = curr->next) {
		currChromSize = hashIntVal(chromSizesHash, curr->chrom);
		if (curr->strand == '+') {
			curr->chromStart = max(0, curr->tss-basalUp);
			curr->chromEnd = min(currChromSize, curr->tss+basalDown);
		}
		else if (curr->strand == '-') {
			curr->chromStart = max(0, curr->tss-basalDown);
			curr->chromEnd = min(currChromSize, curr->tss+basalUp);
		}
		else errAbort("Invalid strand.");
	}

	// Now create the regulatory regions appropriately
	long tmpStart, tmpEnd, prevEnd, nextStart, basalStart, basalEnd;
	for (curr = regdoms; curr != NULL; curr = curr->next) {
		currChromSize = hashIntVal(chromSizesHash, curr->chrom);
		next = curr->next;
		tmpStart = max(0, curr->tss-maximumExtension); // As far as we can extend with no genes nearby
		basalStart = curr->chromStart;
		tmpStart = min(basalStart, tmpStart);
		if ((prev != NULL) && sameString(prev->chrom, curr->chrom)) {
			prevEnd = (prev->strand == '+') ? prev->tss+basalDown : prev->tss+basalUp;
			tmpStart = min(basalStart, max(prevEnd, tmpStart));
		}

		tmpEnd = min(currChromSize, curr->tss + maxExtension);
		basalEnd = curr->chromEnd;
		tmpEnd = max(basalEnd, tmpEnd);
		if ((next != NULL) && sameString(curr->chrom, next->chrom)) {
			nextStart = (next->strand == '+') ? next->tss-basalUp : next->tss-basalDown;
			tmpEnd = max(basalEnd, min(nextStart, tmpEnd));
		}
		curr->chromStart = (int)tmpStart;
		curr->chromEnd = (int)tmpEnd;
		prev = curr;
	}
}

void createOneClosestRegDoms(struct regdom* regdoms, int maximumExtension)
{
	struct regdom *prev = NULL, *curr = NULL, *next = NULL;
	long middle = 0;
	int currChromSize = -1;
	for (curr = regdoms; curr != NULL; curr = curr->next) {
		currChromSize = hashIntVal(chromSizesHash, curr->chrom);
		next = curr->next;
		if ((prev == NULL) || (!sameString(prev->chrom, curr->chrom)))
			prev = NULL;
		if ((next == NULL) || (!sameString(next->chrom, curr->chrom)))
			next = NULL;

		curr->chromStart = max(0, curr->tss - maximumExtension);
		if (prev != NULL) {
			middle = (curr->tss+prev->tss)/2;
			curr->chromStart = max(curr->chromStart, (int)middle);
		}

		curr->chromEnd = min(curr->tss + maximumExtension, currChromSize);
		if (next != NULL) {
			middle = (curr->tss + next->tss)/2;
			curr->chromEnd = min(curr->chromEnd, (int)middle);
		}
		prev = curr;
	}
}

void createTwoClosestRegDoms(struct regdom* regdoms, int maximumExtension)
{
	createBasalPlusExtensionRegDoms(regdoms, maximumExtension, 0, 0);
}


void writeRegulatoryDomains(struct regdom* regdoms, char* fn)
{
	FILE* f = mustOpen(fn, "w");
	struct regdom* curr;
	for (curr = regdoms; curr != NULL; curr = curr->next) {
		// Regulatory domains printed in valid BED6 format, using the following fields
		// chrom    chromStart    chromEnd    name    tss (abuses score field)    strand
		fprintf(f, "%s\t%d\t%d\t%s\t%d\t%c\n", curr->chrom, curr->chromStart, curr->chromEnd, curr->name, curr->tss, curr->strand);
	}
	carefulClose(&f);
}


void createRegulatoryDomains(char* tssFn, char* chromSizesFn, char* association, char* outFn)
/* Create regulatory domains for input genes based on association rule and parameters */
{
	validateInput(association);

	struct regdom* regdoms = readTssFile(tssFn);
	slSort(&regdoms, cmpByChromTssStrand);

	readChromSizes(chromSizesFn);  // Load chrom sizes into global hash

	if (sameString(association, oneClosest))
		createOneClosestRegDoms(regdoms, maxExtension);
	else if (sameString(association, twoClosest))
		createTwoClosestRegDoms(regdoms, maxExtension);
	else if (sameString(association, basalPlusExtension))
		createBasalPlusExtensionRegDoms(regdoms, maxExtension, basalUpstream, basalDownstream);
	else
		errAbort("Association rule must be one of %s, %s, or %s", oneClosest, twoClosest, basalPlusExtension);

	writeRegulatoryDomains(regdoms, outFn);

	regdomFreeList(&regdoms);
	hashFree(&chromSizesHash);
}


int main(int argc, char *argv[])
{
	/* Get all inputs */
	optionInit(&argc, argv, options);
	if (argc != 5) usage();

	char *tssFn, *chromSizesFn, *associationRule, *outFn;
	tssFn = argv[1];
	chromSizesFn = argv[2];
	associationRule = argv[3];
	outFn = argv[4];
	maxExtension = optionInt("maxExtension", maxExtension);
	basalUpstream = optionInt("basalUpstream", basalUpstream);
	basalDownstream = optionInt("basalDownstream", basalDownstream);

	createRegulatoryDomains(tssFn, chromSizesFn, associationRule, outFn);

	optionFree();

	return 0;
}

