/*=========================================================================
 * bdup_mulq.c Compact Givens rotation apply. 
 * bdup_mulq.c calls the Fortran 90 functions bidiag_up_mulq_c,
 * which includes a C interface
 *
 * [X,p] = bdup_mulq(X,Q1,Q2,Q3,trans) computes Givens transformations
 *
 *  X_new(p,:) = Q1*Q2*Q3*X         (trans=1)
 *  X_new(p,:) = Q3'*Q2'*Q1'*X      (trans=0)
 *
 * See bdup_mulp.c for the corresponding applies from the right
 *
 *=========================================================================
 * 01/25/25, J.B., initial interface
 * 03/03/25, J.B., new interface 
 */

#include <string.h> /* needed for memcpy() */
#include "mex.h"
#include "math.h"
#include "bidiag_up_mulq_c.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *X, *Q1, *Q2, *Q3;     /* pointers to input & output matrices*/
    mxInt32 *idx;
#else
    double *X, *Q1, *Q2, *Q3;       /* pointers to input & output matrices*/
    // int *idx;
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

    // Use the transposed method when there are five arguments
    int trans = 0;
    if( nrhs == 5 ){
        trans = 1;
    }

    /* dimensions of input matrices */
    const int m = mxGetM(prhs[0]);  
    const int n = mxGetN(prhs[0]);
    const int n1 = mxGetN(prhs[1]);
    const int n2 = mxGetN(prhs[2]);
    const int n3 = mxGetN(prhs[3]);

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
    plhs[1] = mxCreateNumericMatrix(m, 1, mxINT32_CLASS, mxREAL);
    
#if MX_HAS_INTERLEAVED_COMPLEX
    X       = mxGetDoubles(plhs[0]);
    idx     = mxGetInt32s(plhs[1]);
    
    Q1      = mxGetDoubles(prhs[1]);
    Q2      = mxGetDoubles(prhs[2]);
    Q3      = mxGetDoubles(prhs[3]);
    
#else
    X       = mxGetPr(plhs[0]);
    
    //idx     = mxGetPr(plhs[1]);
    
    int idx_a[m];
    int *idx = idx_a;
    //idx     = mxGetInt32s(plhs[1]);

    Q1      = mxGetPr(prhs[1]);
    Q2      = mxGetPr(prhs[2]);
    Q3      = mxGetPr(prhs[3]);
    
#endif

    // Fortran algorithm

    bidiag_up_mulq_c(X,Q1,Q2,Q3,idx,&m,&n,&n1,&n2,&n3,&trans);

    // Copy values from integer array
#if MX_HAS_INTERLEAVED_COMPLEX
    
#else

    unsigned char *start_of_pr;
    size_t bytes_to_copy;
    start_of_pr = (unsigned char *)mxGetData(plhs[1]);
    bytes_to_copy = m * mxGetElementSize(plhs[1]);
    memcpy(start_of_pr,idx,bytes_to_copy);
    
#endif

    // phase1_c(B1,w,p,Q1,Q2,P1,&m,&n,&m1,&m2);

    // phase2_c(B1,Q3,P2,&m,&n,&m1);
}
