/* Compute the Z-score of the columns of a given matrix.
   If a column has zero standard deviation, its Z-score is considered 0.
19-mar-07	added flag input (ES)
 */

#include <math.h>
#include "mex.h"

/* Input Arguments */
#define	X_IN		prhs[0]

/* Output Arguments */
#define	Z_OUT	0

/*
  X - input matrix of size [n ,p]  (n rows, p columns).
  Z - output matrix of the same size of X.  Memory should have already been
      allocated for Z when compute_fast_zscore() is called. */
void compute_fast_zscore(double *X, unsigned int n, unsigned int p, double f, double *Z)
{
    unsigned int i, j, ind;
    double m, s;
    
    for (i = 0; i < p; i++) {      /* for each column: */

        m = 0;
        /* Compute mean: */
        for (j = 0; j < n; j++) {
            m += *(X + j + i*n);
        }
        m /= n;

        /* compute std. dev.: */
        s = 0;
        for (j = 0; j < n; j++) {
            ind = j + i*n;
            *(Z + ind) = *(X + ind) - m;
            s += (*(Z + ind)) * (*(Z + ind));
        }
        /* s = ((n == 1) ? 0 : sqrt(s / (n - 1))); */
		s = ((n == 1) ? 0 : sqrt(s / (n - f)));

        /* compute Z-score: */
        if (s != 0)
            for (j = 0; j < n; j++) {
                ind = j + i*n;
                *(Z + ind) = *(Z + ind) / s;
            }
        else   /* s == 0 */
            for (j = 0; j < n; j++)
                *(Z + j + i*n) = 0;
            
    }

}

/*
  Usage:
      Z = zsc(X,flag);

  where X is a matrix.
		flag is 1 to normalize by n, 0 to normalize by (n-1)
        Z will be a matrix of the same size of X.

  Z-score is computed for each columns separately.
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* dimensions. */
    double *X;
    double *Z;
    unsigned int n, p;
	double f;

    /* Check for proper number of arguments */
    if (nrhs != 2)
	mexErrMsgTxt("zsc requires two input arguments.");
    else if (nlhs > 1)
	mexErrMsgTxt("zsc requires one output arguments.");
    
    /* Get input dimensions. */
    n = mxGetM(X_IN);
    if (n == 0)
	mexErrMsgTxt("X should have at least 1 row.");
    p = mxGetN(X_IN);
    if (p == 0)
	mexErrMsgTxt("X should have at least 1 column.");

    /* Assign pointers to the various parameters */
    /* /X = mxGetPr(X_IN); */
	X = mxGetPr(prhs[0]);	/* / data */
	/* /f = mxGetPr(prhs[1]);	// flag */
	f = mxGetScalar(prhs[1]);
	f = ( f == 0 ) ? 1 : 0;	/* / if zero, normalize by n-1 */

    /* Allocate memory and assign output pointer */
    plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL); /* mxReal is our data-type */
    Z = mxGetPr(plhs[0]);

    /* Do the actual computations in a subroutine */
    compute_fast_zscore(X, n, p, f, Z);
}

