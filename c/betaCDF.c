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

/* constants for Lentz's method for computing the incomplete beta function
 */
#define MAXIT 10000
#define EPS 3.0e-7 
#define FPMIN 1.0e-30 

struct optionSpec options[] = {
	{NULL, 0}
};

void usage()
{
  errAbort(
  "\n"
  "Usage: betaCDF alpha beta x\n"
  );
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

int main(int argc, char *argv[])
{
	/* Get all inputs */
	optionInit(&argc, argv, options);
	if (argc != 4) usage();

	double alpha, beta, x;
	alpha = atof(argv[1]);
	beta = atof(argv[2]);
	x = atof(argv[3]);

	printf("%f\n", betai(alpha, beta, x));

	optionFree();

	return 0;
}
