/* /////////////////////////////// CALC_COR.C /////////////////
//
// 01-nov-03 ES */

#include "mex.h"

/* calculate delays by enumerating on first spike train only */

void calc_delays(double *s1, int len1, double *s2, int len2, int OS, double *h)
{
	int d, d1, len, p1, p2, p20, count;

	len = 2 * OS + 1;
	p1 = 0;
	p2 = 0;
	p20 = 0;
	count = 0;

	/* printf("CALC_DC.C init: len1 = %d; len2 = %d; len = %d\n", len1, len2, len ); */

	while ( p1 < len1 ) {

		d = (int)(s2[ p2 ] - s1[ p1 ]);
		if ( d < (-OS) ) {						/* s2 too early */
			p20++;
			if ( p20 == len2 ) break;
			p2 = p20;
		}
		else {
			if ( d > OS ) {						/* s2 too late */
				p1++;
				if ( p1 == len1 ) break;
				d1 = (int)(s1[ p1 ] - s1[ p1 - 1 ]);
				if ( d1 < len )
					p2 = p20;
			}
			else {								/* in range */
				(h[ d + OS ])++;
				p2++;
				if ( p2 == len2 ) {
					p2 = p20;
					p1++;
				}
			}
		}
		count++;

	}

	/* printf("CALC_DC.C count: %d; p1 = %d; p2 = %d\n", count, p1, p2 ); */

}

/* gateway */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
  double offset;
  double *vec1, *vec2, *corrVec;
  int lvec1,lvec2,m,n;
  
  /* check io */
  if (nrhs != 3) mexErrMsgTxt("exactly 3 arguments required (sp1, sp2, offset)");
  if (nlhs > 1) mexErrMsgTxt("only 1 output required (corr)");

  /* get input */
  vec1 = mxGetPr(prhs[0]);
  lvec1 = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]); 
  if ( (n>1) || mxIsSparse(prhs[0]) )
	  mexErrMsgTxt("sp1 must be a full column vector");

  vec2 = mxGetPr(prhs[1]);
  lvec2 = mxGetM(prhs[1]);
  n = mxGetN(prhs[1]); 
  if ( (n>1) || mxIsSparse(prhs[1]) )
	  mexErrMsgTxt("sp2 must be a full column vector");

  m = mxGetM(prhs[2]);
  n = mxGetN(prhs[2]);
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(m == 1 && n == 1)) {
    mexErrMsgTxt("Offset must be a noncomplex scalar double.");
  }
  offset = mxGetScalar(prhs[2]);
  if (offset<0) mexErrMsgTxt("offset must be non-negative");

  /* allocate output */
  plhs[0] = mxCreateDoubleMatrix(offset*2+1, 1, mxREAL);
  corrVec = mxGetPr(plhs[0]);
  
  /* compute */
  calc_delays(vec1, lvec1, vec2, lvec2, offset, corrVec);
}

