/*=========================================================
 * bidiagup.c Bidiagonal updating algorithm. 
 * bidiagup.c calls the Fortran 90 functions phase1_c and phase2_c,
 * which include a C interface
 *
 * [B1,Q1,Q2,Q3,P1,P2] = bidiag(B,w,p) computes the bidiagonal in B1 and
 *  Q1,Q2,Q3 and P1,P2 which compactly store the givens rotations 
 *
 * This is a MEX-file for MATLAB.
 *
 *=========================================================*
 * 01/25/25, J.B., initial interface
 */

#include "mex.h"
#include "math.h"
#include "phase1_c.h"
#include "phase2_c.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *B1, *w, *p, *Q1, *Q2, *Q3, *P1, *P2; /* pointers to input & output matrices*/
#else
    double *B1, *w, *p, *Q1, *Q2, *Q3, *P1, *P2; /* pointers to input & output matrices*/
#endif
    // ptrdiff_t m,n;      /* matrix dimensions */
    // int m1;
    
    /*
#if MX_HAS_INTERLEAVED_COMPLEX
    B = mxGetDoubles(prhs[0]); 
    w = mxGetDoubles(prhs[1]); 
    p = mxGetDoubles(prhs[2]); 
#else
    B = mxGetPr(prhs[0]); 
    w = mxGetPr(prhs[1]); 
    p = mxGetPr(prhs[2]); 
#endif
*/
    /* dimensions of input matrices */
    const int m = mxGetM(prhs[0]);  
    const int n = mxGetN(prhs[0]);
    const int m1 = ceil((n*n)/2.0)+1;   /* max number of rotations */
    const int m2 = (m-n)+1;

    /*
    // check to make sure the first input argument is a real matrix
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
      mexErrMsgIdAndTxt( "MATLAB:matrixMultiply:fieldNotRealMatrix",
              "First input argument must be a real matrix.");
    }

    // check to make sure the second input argument is a real matrix //
    if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
      mexErrMsgIdAndTxt( "MATLAB:matrixMultiply:fieldNotRealMatrix",
              "Second input argument must be a real matrix.");
    }

    if (p != mxGetM(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
                "Inner dimensions of matrix multiply do not match.");
    }
    */

    /* create output matrices B4, Q, Q1, P and also inputs */
    plhs[0] = mxDuplicateArray(prhs[0]);
    // plhs[1] = mxCreateDoubleMatrix(m1, 4, mxREAL);
    // plhs[2] = mxCreateDoubleMatrix(m2, 4, mxREAL);
    // plhs[3] = mxCreateDoubleMatrix(m1, 4, mxREAL);
    // plhs[4] = mxCreateDoubleMatrix(m1, 4, mxREAL);
    // plhs[5] = mxCreateDoubleMatrix(m1, 4, mxREAL);

    plhs[1] = mxCreateDoubleMatrix(4, m1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(4, m2, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(4, m1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(4, m1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(4, m1, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    B1      = mxGetDoubles(plhs[0]);
    Q1      = mxGetDoubles(plhs[1]);
    Q2      = mxGetDoubles(plhs[2]);
    Q3      = mxGetDoubles(plhs[3]);
    P1      = mxGetDoubles(plhs[4]);
    P2      = mxGetDoubles(plhs[5]);

    w       = mxGetDoubles(mxDuplicateArray(prhs[1]));
    p       = mxGetDoubles(mxDuplicateArray(prhs[2]));
#else
    B1      = mxGetPr(plhs[0]);
    Q1      = mxGetPr(plhs[1]);
    Q2      = mxGetPr(plhs[2]);
    Q3      = mxGetPr(plhs[3]);
    P1      = mxGetPr(plhs[4]);
    P2      = mxGetPr(plhs[5]);

    w       = mxGetPr(mxDuplicateArray(prhs[1]));
    p       = mxGetPr(mxDuplicateArray(prhs[2]));
#endif

    // Fortran algorithms

    phase1_c(B1,w,p,Q1,Q2,P1,&m,&n,&m1,&m2);

    phase2_c(B1,Q3,P2,&m,&n,&m1);
}
