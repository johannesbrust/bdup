/*=========================================================
 * bhu.c Bidiagonal updating algorithm. 
 * bhu.c calls the Fortran 90 functions bhu_c, which is a 
 * C interface to a "sparse" householder algorithm
 *
 * [B1,Y,W] = bhu(B,w,p,kk) computes the bidiagonal in B1 where
 * Y and W compactly store the Householder reflectors. The
 * method can compute a truncated factorization, for kk
 * iterations.
 *
 * The input B and output B1 is a set of two columns
 *
 * This is a MEX-file for MATLAB.
 *
 *=========================================================*
 * 01/25/25, J.B., initial interface
 */

#include "mex.h"
#include "math.h"
#include "bhu_c.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *Y, *W, *B, *BO, *wh, *p; /* pointers to input & output matrices*/
#else
    double *Y, *W, *B, *BO, *wh, *p; /* pointers to input & output matrices*/
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
    const int m = mxGetM(prhs[1]);  
    const int n = mxGetM(prhs[2]);
    //const int kk = prhs[3];

    const int kk = (int)mxGetScalar(prhs[3]);
    
    //const int m1 = ceil((n*n)/2.0)+1;   /* max number of rotations */
    //const int m2 = (m-n)+1;

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

    /* create output matrices BO, Y, and W*/
    plhs[0] = mxCreateDoubleMatrix(kk, 2, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(m, kk, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, kk-1, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    BO      = mxGetDoubles(plhs[0]);
    Y       = mxGetDoubles(plhs[1]);
    W       = mxGetDoubles(plhs[2]);

    B       = mxGetDoubles(prhs[0]);
    wh      = mxGetDoubles(prhs[1]);
    p       = mxGetDoubles(prhs[2]);
#else
    BO      = mxGetPr(plhs[0]);
    Y       = mxGetPr(plhs[1]);
    W       = mxGetPr(plhs[2]);
    
    B       = mxGetPr(prhs[0]);
    wh      = mxGetPr(prhs[1]);
    p       = mxGetPr(prhs[2]);
#endif

    // Fortran algorithms

    bhu_c(Y,W,B,BO,wh,p,&m,&n,&kk);

}
