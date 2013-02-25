#include "regdom.h"
#include "common.h"
#include "linefile.h"

int cmpByChromTssStrand(const void* pa, const void* pb)
/* Compare two regdoms based on chromosome, tss, and strand */
{
	const struct regdom* a = *((struct regdom**)pa);
	const struct regdom* b = *((struct regdom**)pb);
	int chromCmp = strcmp(a->chrom, b->chrom);
	if (chromCmp < 0) return -1;
	else if (chromCmp > 0) return 1;

	if (a->tss < b->tss) return -1;
	else if (a->tss > b->tss) return 1;

	if ((a->strand == '+') && (b->strand == '-')) return -1;
	else if ((a->strand == '-') && (b->strand == '+')) return 1;

	return strcmp(a->name, b->name);
}

void regdomFree(struct regdom** pEl)
/* Free a dynamically allocated regdom */
{
	struct regdom* el;
	if ((el = *pEl) == NULL) return;
	freeMem(el->chrom);
	freeMem(el->name);
	freez(pEl);
}

void regdomFreeList(struct regdom** pList)
/* Free a list of dynamically allocated regdoms */
{
	struct regdom *el, *next;
	for (el = *pList; el != NULL; el = next) {
		next = el->next;
		regdomFree(&el);
	}
	*pList = NULL;
}

struct regdom* readTssFile(char* fn)
/* Reads a file of "chrom \t tss \t strand \t name" elements into a list of regulatory domains */
{
	int fieldCount = 4;
	char* row[fieldCount];
	struct regdom* regdomList = NULL;
	int wordCount;
	struct lineFile* lf = lineFileOpen(fn, TRUE);
	while ((wordCount = lineFileChopTab(lf, row)) != 0) {
		if (wordCount != fieldCount)
			errAbort("Expecting exactly %d words line %d of %s", fieldCount, lf->lineIx, lf->fileName);
		struct regdom* el;
		AllocVar(el);
		el->chrom = cloneString(row[0]);
		el->tss = lineFileNeedNum(lf, row, 1);
		el->strand = row[2][0];
		el->name = cloneString(row[3]);
		el->chromStart = el->tss;
		el->chromEnd = el->tss;
		slAddHead(&regdomList, el);
	}
	lineFileClose(&lf);
	return regdomList;
}

struct regdom* readInitializedRegdomFile(char* fn)
/* Reads a file of "chrom \t chromStart \t chromEnd \t name \t tss \t strand" elements (initialized regdoms) into
 a linked list of regulatory domains */
{
	int fieldCount = 6;
	char* row[fieldCount];
	struct regdom* regdomList = NULL;
	int wordCount;
	struct lineFile* lf = lineFileOpen(fn, TRUE);
	while ((wordCount = lineFileChopTab(lf, row)) != 0) {
		if (wordCount != fieldCount)
			errAbort("Expecting exactly %d words line %d of %s", fieldCount, lf->lineIx, lf->fileName);
		struct regdom* el;
		AllocVar(el);
		el->chrom = cloneString(row[0]);
		el->chromStart = lineFileNeedNum(lf, row, 1);
		el->chromEnd = lineFileNeedNum(lf, row, 2);
		el->name = cloneString(row[3]);
		el->tss = lineFileNeedNum(lf, row, 4);
		el->strand = row[5][0];
		slAddHead(&regdomList, el);
	}
	lineFileClose(&lf);
	return regdomList;
}

