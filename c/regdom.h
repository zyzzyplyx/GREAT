#ifndef REGDOM_H
#define REGDOM_H

struct regdom {
	struct regdom* next;
	char* chrom;
	int chromStart;     // The leftmost edge of the regulatory domain
	int chromEnd;       // The rightmost edge of the regulatory domain
	int tss;
	char strand;
	char* name;         // The unique identifier of the gene
};

/* Comparison function */
int cmpByChromTssStrand(const void* pa, const void* pb);

/* Free a dynamically allocated regdom */
void regdomFree(struct regdom** pEl);

/* Free a list of dynamically allocated regdoms */
void regdomFreeList(struct regdom** pList);

/* Reads a file of "chrom \t tss \t strand \t name" elements (uninitialized regdoms) into a list of regulatory domains */
struct regdom* readTssFile(char* fn);

/* Reads a file of "chrom \t chromStart \t chromEnd \t name \t tss \t strand" elements (initialized regdoms) into
 a linked list of regulatory domains */
struct regdom* readInitializedRegdomFile(char* fn);

#endif
