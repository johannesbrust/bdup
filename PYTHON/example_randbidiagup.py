##########################################################################################
#
# Test computing the bidiagonal factorization of an updated matrix
# 
########################################################################################## 
# 02/21/25, J.B., initial version

import numpy as np
from bgu import bgu
import time

def main():
    
    m = 50000
    n = 3000
    seed_bd = 1

    # Initialize random number generator
    
    td  = time.time() 
    rng = np.random.default_rng(seed_bd)

    B   = np.diag(rng.standard_normal(n) ) + np.diag(rng.standard_normal(n-1),1) #rng.standard_normal((m,n)) * np.tri()
    B   = np.vstack((B,np.zeros((m-n,n))))
    w   = rng.standard_normal(m)
    p   = rng.standard_normal(n)
    wcp = w.copy()
    pcp = p.copy()

    tde = time.time() - td

    # Bidiagonal updating algorithm
    #

    ta = time.time()

    B1, Q1, Q2, Q3, P1, P2 = bgu(B,w,p)

    tae = time.time() - ta

    err = np.abs(np.linalg.norm(B1,'fro') - np.linalg.norm(B+np.linalg.outer(wcp,pcp),'fro'))

    print('');
    print('************** Bidiagonal Update ********************  ');
    print('*         Alg: bgu (Python 3.9)                        ');
    print('*         Size:      m = %i, n = %i                    '% (m,n));
    print('*         Error      = %1.8f                           '% err);
    print('*         Time data  = %1.5f                           '% tde);
    print('*         Time alg   = %1.5f                           '% tae);        
    print('*****************************************************  ');        
    print('');

if __name__ == "__main__":
    main()