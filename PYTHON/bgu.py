##########################################################################################
#
# Interface to fortran bidiagonal updating algorithms
#
#   B1, Q1, Q2, Q3, P1, P2 = bidiagup(B,w,p)
#
# Computes an updated bidiagonal matrix from a previous bidiagonal B, and rank-1 update
# w * p' in the new B1 and orthognal rotation information in Q1, Q2, Q3 and P1 and P2
# If Qi, and Pi were assembled, which they are not, the computation corresponds to
#
#   B1 = Q3'*Q2'*Q1'*(B+w*p')*P1*P2
# 
# This method is for m >= n where m=nrows(B) and n = ncols(B)
#
# Version Aug. 29, 2025
# J.J. Brust, 
# johannesbrust@yahoo.com
##########################################################################################
# 02/19/25, J.B., Initial version
# 02/21/25, J.B., Linux interface
# 02/24/25, J.B., Developing the interface
# 03/05/25, J.B., Copies of low rank update vectors, and givens multiply functions
# 08/29/25, J.B., Conditional loading of library
# 08/30/25, J.B., Improving copy operations and prepare for release
# 08/31/25, J.B., Refactor and prepare for release

import ctypes
import numpy as np
import platform
from ctypes import byref, c_int
from numpy.ctypeslib import load_library, as_array
from pathlib import Path

# Load library
lib_path1 = Path(__file__).parent.parent / "FORTRAN" / "objs"

if platform.machine() == 'x86_64':
    if platform.system() == 'Darwin':
        lib_path = lib_path1 / "osx_x64"
    elif platform.system() == 'Linux':
        lib_path = lib_path1 / "unix_x64"
elif platform.machine() == 'arm64':
    lib_path = lib_path1 / "osx_arm"

#print(str(lib_path))

lib = load_library("bidiagup", str(lib_path))

def bgu(B,w,p,dat='C'):

    # Python memory is row-wise, but Fortran has a column-wise layout
 
    # Dimensions
    m = np.size(w)
    n = np.size(p)

    # Transpose the data if dimension does not correspond to m > n
    if n > m: 
        w_tmp = w
        w = p
        p = w_tmp
        B = B.T
        m = np.size(w)
        n = np.size(p)

    # Avoiding the copy operation by constructing a new bidiagonal from the
    # old diagonal and upper diagonal
    B1 = np.zeros((m,n),order='F')
    np.fill_diagonal(B1,np.diagonal(B).copy())
    np.fill_diagonal(B1[:-1, 1:], np.diagonal(B,1).copy())
    
    # Copy w and p (two vectors are assumped to be relatively inexpensive)
    wcp     = w.copy()
    pcp     = p.copy()
    
    m1 = (np.ceil(n*n/2.0)+1).astype(np.int64)
    m2 = m-n+1

    # Compact reflector data
    Q1      = np.zeros((4,m1),order='F') 
    Q2      = np.zeros((4,m2),order='F')
    Q3      = np.zeros((4,m1),order='F')
    P1      = np.zeros((4,m1),order='F')
    P2      = np.zeros((4,m1),order='F')

    # C types    
    c_B     = cpntr(B1)
    c_w     = cpntr(wcp)
    c_p     = cpntr(pcp)       
    c_Q1    = cpntr(Q1)
    c_Q2    = cpntr(Q2)
    c_Q3    = cpntr(Q3)
    c_P1    = cpntr(P1)
    c_P2    = cpntr(P2)
    c_m     = c_int(m)
    c_n     = c_int(n)
    c_m1    = c_int(m1)
    c_m2    = c_int(m2)

    #
    # Call algorithm
    # 
    lib.phase1_c(c_B, c_w, c_p, c_Q1, c_Q2, c_P1, byref(c_m), byref(c_n), byref(c_m1), byref(c_m2))
    lib.phase2_c(c_B, c_Q3, c_P2, byref(c_m), byref(c_n), byref(c_m1), byref(c_m2))

    # Return values
    B1 = as_array(c_B,  shape=(m*n,) ).reshape((m,n),  order='F')    
    Q1 = as_array(c_Q1, shape=(4*m1,)).reshape((4,m1), order='F')
    Q2 = as_array(c_Q2, shape=(4*m2,)).reshape((4,m2), order='F')
    Q3 = as_array(c_Q3, shape=(4*m1,)).reshape((4,m1), order='F')
    P1 = as_array(c_P1, shape=(4*m1,)).reshape((4,m1), order='F')
    P2 = as_array(c_P2, shape=(4*m1,)).reshape((4,m1), order='F')

    if dat == 'C':
        return B1.copy(order='C'), Q1.copy(order='C'), Q2.copy(order='C'), Q3.copy(order='C'), \
               P1.copy(order='C'), P2.copy(order='C')
    else:
        return B1, Q1, Q2, Q3, P1, P2

def bgu_mulq(X,Q1,Q2,Q3,trans,dat='C'):

    # Python memory is row-wise, but Fortran has a column-wise layout

    # Dimensions
    m,n = X.shape

    n1  = Q1.shape[1] 
    n2  = Q2.shape[1]
    n3  = Q3.shape[1]
        
    X1  = np.array(X, order='F', copy=True)
    idx = np.asarray(np.arange(0,m))
    
    # C types
    c_X     = X1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_idx   = idx.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    
    # Pointers
    # c_Q1    = ens_fort(Q1).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # c_Q2    = ens_fort(Q2).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # c_Q3    = ens_fort(Q3).ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    c_Q1    = cpntr(Q1)
    c_Q2    = cpntr(Q2)
    c_Q3    = cpntr(Q3)
    
    c_m     = c_int(m)
    c_n     = c_int(n)
    c_n1    = c_int(n1)
    c_n2    = c_int(n2)
    c_n3    = c_int(n3)

    c_trans = c_int(trans)

    #
    # Call algorithm
    # 
    lib.bidiag_up_mulq_c(c_X,c_Q1,c_Q2,c_Q3,c_idx,byref(c_m),byref(c_n),byref(c_n1),byref(c_n2),byref(c_n3),byref(c_trans));

    # Return values

    X1  = as_array(c_X, shape=(m*n,)).reshape((m,n), order='F')
    idx = as_array(c_idx,idx.shape) -1

    if dat == 'C':
        return X1.copy(order='C'), idx 
    else:
        return X1, idx
        
def bgu_mulp(X,P1,P2,trans,dat='C'):

    # Python memory is row-wise, but Fortran has a column-wise layout

    # Dimensions
    m,n = X.shape

    n1  = P1.shape[1] 
    n2  = P2.shape[1]    
    
    X1 = np.array(X, order='F', copy=True)  
    idx = np.asarray(np.arange(0,n))
    
    # C types

    c_X     = X1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_idx   = idx.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    
    # Pointers
    
    # c_P1    = ens_fort(P1).ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    # c_P2    = ens_fort(P2).ctypes.data_as(ctypes.POINTER(ctypes.c_double))    
    
    c_P1    = cpntr(P1)
    c_P2    = cpntr(P2)

    c_m     = c_int(m)
    c_n     = c_int(n)
    c_n1    = c_int(n1)
    c_n2    = c_int(n2)    

    c_trans = c_int(trans)

    #
    # Call algorithm (Fortran algorithm via C interface)
    # 
    lib.bidiag_up_mulp_c(c_X,c_P1,c_P2,c_idx,byref(c_m),byref(c_n),byref(c_n1),byref(c_n2),byref(c_trans));

    X1  = as_array(c_X, shape=(m*n,)).reshape((m,n), order='F')
    idx = as_array(c_idx,idx.shape) -1

    if dat == 'C':
        return X1.copy(order='C'), idx
    else:
        return X1, idx 
    
#
# Auxiliary functions
#

# ENS_FORT(A) Ensures that A is in 
# Fortran column major memory 
def ens_fort(A):
    if not A.flags['F_CONTIGUOUS']:
        return np.array(A, order='F', copy=True)
    return A

# CPNTR(A) Creates a C-Pointer to 
# Fortran contiguous data 
def cpntr(A):
    return ens_fort(A).ctypes.data_as(ctypes.POINTER(ctypes.c_double))