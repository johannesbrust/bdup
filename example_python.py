# example.py
"""
Example of the bidiagonal updating algorithms: bgu and bhu

This script computes the bidiagonal factorization of B + w * p.T.
Depending on the algorithm, the orthogonal matrices may be formed explicitly.

Author: J.J. Brust (johannesbrust@yahoo.com)
Date: Initial 04/08/25, Prep for release 08/26/25, Final prep 08/27/25
"""

import numpy as np
import time
import platform
from PYTHON.bgu    import bgu, bgu_mulq, bgu_mulp
from PYTHON.bhu    import bhu

if __name__ == "__main__":

    # Setup problem
    fc = 200        # 1
    m = fc * 10
    n = fc * 4
    B = np.vstack([
        np.diag(np.arange(1, n+1)) + np.diag(np.arange(2, n+1), 1),
        np.zeros((m - n, n))
    ])

    w = np.zeros((m, 1))
    p = np.zeros((n, 1))
    w[1] = 1
    w[5] = 10
    p[1] = 2
    p[2] = 3

    if m >= 500:
        w[499] = 1

    print("                                                   ")
    print("***************** Bidiag update *******************")
    print("*         Algorithms: bgu                          ")
    print("*                     bhu                          ")
    print("*                                                  ")
    print(f"*           Software: Python {platform.python_version()} ")
    print(f"*             System: {platform.machine()}        ")
    print("*                                                  ")
    print(f"*          Prob size: m={m}")
    print(f"*                     n={n}")
    print(f"*                     nnz(w)={np.count_nonzero(w)}")
    print(f"*                     nnz(p)={np.count_nonzero(p)}")
    print("*                                                  ")
    print("*         Release Aug. 2025                        ")
    print("*         J.J. Brust (johannesbrust@yahoo.com)     ")
    print("***************************************************\n")

    # BGU algorithm (Givens-based)
    start_bgu = time.time()
    B1, Q1, Q2, Q3, P1, P2 = bgu(B, w, p, 'F')
    tbgu = time.time() - start_bgu

    # Error or full reconstruction
    if m <= 1000:
        start_mult_bgu = time.time()
        BU, pq = bgu_mulq(B1, Q1, Q2, Q3, 1)
        BU, pp = bgu_mulp(BU, P1, P2, 1)
        tmbgu = time.time() - start_mult_bgu

        B_orig = B + w @ p.T
        err = np.linalg.norm(BU[pq][:, pp] - B_orig, ord='fro')
    else:
        tmbgu = tbgu
        err = abs(np.linalg.norm(B1, ord='fro') - np.linalg.norm(B + w @ p.T, ord='fro'))

    # BHU algorithm (Householder-based)
    start_bhu = time.time()
    BUP = np.column_stack([
        np.diag(B[:n, :n]),
        np.append(np.diag(B[:n, :n], k=1), 0)
    ])
    B1h, Y, W = bhu(BUP, w, p, n)
    tbhu = time.time() - start_bhu

    if m <= 1000:
        start_mult_bhu = time.time()
        YL = np.tril(Y)
        WL = np.tril(W)

        cls_Y = Y.shape[1]
        cls_W = W.shape[1]

        T = np.triu(Y[:cls_Y, :cls_Y], 1) + np.eye(cls_Y)
        R = np.triu(W[:cls_W, :cls_W], 1) + np.eye(cls_W)

        QQ = np.eye(m) - 2 * YL @ np.linalg.solve(T, YL.T)
        PP = np.eye(n) - 2 * WL @ np.linalg.solve(R, WL.T)

        B1h_diag = np.diag(B1h[:n, 0])
        B1h_super = np.diag(B1h[:n-1, 1], k=1)
        B1h_full = B1h_diag + B1h_super

        errh = np.linalg.norm(B + w @ p.T - QQ[:, :n] @ B1h_full @ PP[:, :n].T, ord='fro')
        tmbhu = time.time() - start_mult_bhu
    else:
        tmbhu = tbhu
        errh = abs(np.linalg.norm(B1h, ord='fro') - np.linalg.norm(B + w @ p.T, ord='fro'))

    # Results
    print(" alg: bgu")
    print(f" time update: {tbgu:.4f} ")
    print(f" time mult:   {tmbgu:.4f}")
    print(f" error:       {err:.5e}  ")
    print()
    print(" alg: bhu")
    print(f" time update: {tbhu:.4f} ")
    print(f" time mult:   {tmbhu:.4f}")
    print(f" error:       {errh:.5e} ")
    print("                          ")
