##########################################################################################
#
# Interface to Fortran bidiagonal updating algorithms
#
#   B1, Y, W = bhu(B,w,p,kk,dat)
#
# Computes the bidiagonal in B1 where Y and W compactly store the Householder reflectors. 
# The method can compute a truncated factorization, for kk iterations. If dat='C',
# the default, the result is in row-wise memory layout. Setting dat='F' is faster,
# and returns the arrays in column majour layout.
#
# The input B and output B1 is a set of two columns# 
# This method is for m >= n where m=nrows(B) and n = ncols(B)
#
# Version Aug. 29, 2025
# J.J. Brust, 
# johannesbrust@yahoo.com
##########################################################################################
# 08/29/25, J.B., Initial version of bhu interface
# 08/30/25, J.B., Optionally return data in column-major order (this is faster)

import ctypes
import numpy as np
import platform
from ctypes import byref, c_int
from numpy.ctypeslib import load_library, as_array
from pathlib import Path

# Set trunk of path
lib_path1 = Path(__file__).parent.parent / "FORTRAN" / "objs"

# Load library depending on architecture and OS
if platform.machine() == 'x86_64':
    if platform.system() == 'Darwin':
        lib_path = lib_path1 / "osx_x64"
    elif platform.system() == 'Linux':
        lib_path = lib_path1 / "unix_x64"
elif platform.machine() == 'arm64':
    lib_path = lib_path1 / "osx_arm"

#print(str(lib_path))

lib = load_library("bhu", str(lib_path))

def bhu(B,w,p,kk=float('inf'),dat='C'):

    # Python memory is row-wise, but Fortran has a column-wise layout
    # This incurrs potentially expensive copy operations
    
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

    # Possibly truncate factorization 
    kk = int(np.minimum(n,kk))

    # B is two columns (copy is expected to be rel. inexpensive)    
    Bcp = B.copy(order='F')
    B1  = np.zeros((kk,2), order='F')

     # Copy w and p (two vectors are assumped to be relatively inexpensive)
    wcp     = w.copy()
    pcp     = p.copy()
    
    # Initialize W and Y
    Y = np.zeros( (m,kk), order='F' )
    W = np.zeros( (n,kk), order='F' )

    # C types
    c_B     = B1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_Bcp   = Bcp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_w     = wcp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_p     = pcp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_Y     = Y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    c_W     = W.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
 
    c_m     = c_int(m)
    c_n     = c_int(n)
    c_kk    = c_int(kk)
    
    #
    # Call algorithm
    # 
    lib.bhu_c(c_Y,c_W,c_Bcp,c_B,c_w,c_p,byref(c_m),byref(c_n),byref(c_kk))

    # Form numpy arrays
    Y  = as_array(c_Y, shape=(kk * m,)).reshape((m, kk), order='F')
    W  = as_array(c_W, shape=(kk * n,)).reshape((n, kk), order='F')
    B1 = as_array(c_B, shape=(kk * 2,)).reshape((kk, 2), order='F')

    if dat == 'C':
        return B1, Y, W
    else:
        return B1.copy(order='C'), Y.copy(order='C'), W.copy(order='C')

    # Row-major order
    # B1 = B1_f.copy(order='C')
    # Y  = Y_f.copy(order='C')
    # W  = W_f.copy(order='C')
    
    # return B1, Y, W
