##########################################################################################
#
# Test for updating a svd using the bidiagonal updating algorithm and the build-in svd
# functions in SciPy
# 
########################################################################################## 
# 02/21/25, J.B., initial version
# 03/05/25, J.B., svd version
# 08/30/25, J.B., preparation for release (incl. profiling)

import numpy as np
#from bidiagup import bidiagup
from bgu import bgu
import time
import scipy as sc

import cProfile
import pstats

def main():
    
    fac = 1             # 1
    fac1= 10000         # 1
    fac2= 1000          # 1
    m   = fac*fac1      # 10000
    n   = fac*fac2      # 1000
    seed_bd = 1

    # Initialize random number generator
    
    td  = time.time() 
    rng = np.random.default_rng(seed_bd)

    #B   = np.diag(rng.standard_normal(n) ) + np.diag(rng.standard_normal(n-1),1) #rng.standard_normal((m,n)) * np.tri()
    #   = np.vstack((B,np.zeros((m-n,n))))
    
    A   = rng.standard_normal((m,n))    
    w   = rng.standard_normal(m)
    p   = rng.standard_normal(n)
    wcp = w.copy()
    pcp = p.copy()

    tde = time.time() - td

    # Initial SVD computation

    tv  = time.time()
    B   = sc.linalg.svdvals(A) 
    Bcp = np.vstack((np.diag(B.copy()),np.zeros((m-n,n))))
    Bcp1= Bcp.copy()
    tve = time.time() - tv

    print('');    
    print('*         Alg: svdvals (scipy)                       ')
    print('*         m = %i, n = %i                             '% (m,n));
    print('*         Time alg. = %1.8f                          '% tve);
    print(''); 

    #
    # Bidiagonal updating algorithm
    #

    ta = time.time()

    profiler = cProfile.Profile()
    profiler.enable()
    
    B1, Q1, Q2, Q3, P1, P2 = bgu(Bcp,w,p,'F')

    profiler.disable()

    tae = time.time() - ta

    # Form a tridiagonal matrix
    tms  = time.time()
    b1   = np.diag(B1)
    b2   = np.diag(B1,1)
    t1   = b1 * b1
    t1[1:] = t1[1:] + b2[:]*b2[:]
    t2 = b1[:n-1]*b2[:]

    #T    = np.matmul(np.transpose(B1),B1)
    tme  = time.time() - tms

    # DEBUG 
    #print(t1)
    #print(t2)
    #print(np.matmul(np.transpose(B1),B1))

    # Compute eigenvalues of the tridiagonal
    tes   = time.time()
    #evals = sc.linalg.eigh_tridiagonal( np.diag(T), np.diag(T,1), eigvals_only=True )
    evals = sc.linalg.eigh_tridiagonal( t1, t2, eigvals_only=True )
    tee   = time.time() - tes
    
    #
    # Direct algorithm
    #
    tss   = time.time()
    svdd  = sc.linalg.svdvals( Bcp1 + np.outer(wcp,pcp) )
    tse   = time.time() - tss

    #err = np.abs(np.linalg.norm(B1,'fro') - np.linalg.norm(B+np.linalg.outer(wcp,pcp),'fro'))

    err = np.linalg.norm(svdd - np.sqrt(evals[::-1]))

    print('');
    print('********************* SVD Update ********************* ');
    print('*                                                      ');
    print('*  Algs: bgu (Python 3.9)                              ');
    print('*        eigh_tridiagonal (SciPy 1.15.1)               ');
    print('*        svdvals (SciPy 1.15.1)                        ');
    print('*                                                      ');
    print('*        Size:              m = %i, n = %i             '% (m,n));
    print('*        Error              = %1.8f                    '% err);
    print('*                                                      ');
    print('*  Time: data               = %1.5f                    '% tde);
    print('*                                                      ');
    print('*  Time: initial svd        = %1.5f                    '% tve);
    print('*                                                      ');
    print('*  Time: update alg                                    ');
    print('*        bgu                = %1.5f                    '% tae);
    print('*        matmul             = %1.5f                    '% tme);
    print('*        eigh_tridiagonal   = %1.5f                    '% tee);
    print('*        total (update alg) = %1.5f                    '% (tae+tme+tee));
    print('*                                                      ');
    print('*  Time: direct alg         = %1.5f                    '% tse);        
    print('*                                                      ');
    print('****************************************************** ');        
    print('');

    # DEBUG
    #print(svdd)
    #print(np.sqrt(evals))

    #print(Bcp)
    #print(wcp)
    #print(pcp)

    # Profiler output

    stats = pstats.Stats(profiler)
    stats.sort_stats('cumtime') # Sort by cumulative time
    stats.print_stats()

if __name__ == "__main__":
    main()