## Bidiagonal updating algorithms

Installation instructions for two bidiagonal updating algorithms both
located in FORTRAN/include/ `[1]`

    bidiag_up.f90           (givens)
    bidiagh_up.f90          (compact householder)

The first algorithm uses Givens rotations to restore the updated bidiagonal
to new bidiagonal again. The second algorithm uses a compact representation of
Householder reflectors.

The algorithms are typically "thin", i.e. economic, so that
orthogonal matrices are not explicitly computed but represented in factored
forms. MATLAB and Python interfaces are provided.

We ship precompiled binaries for linux 64bit x86, mac 64bit x86 and mac 64bit arm.
These have been tested on Apple Silicon (Apple M2 Max) with Python 3.10.18 and 
Matlab R2023b and Intel chips with Python 3 and Matlab . Once the installation is 
complete you can run an example in Python like:

``python python_example.py``

```
***************** Bidiag update *******************
*         Algorithms: bgu                          
*                     bhu                          
*                                                  
*           Software: Python 3.10.18 
*             System: x86_64        
*                                                  
*          Prob size: m=2000
*                     n=800
*                     nnz(w)=3
*                     nnz(p)=2
*                                                  
*         Release Aug. 2025                        
*         J.J. Brust (johannesbrust@yahoo.com)     
***************************************************

 alg: bgu
 time update: 0.0671 
 time mult:   0.0671
 error:       3.63798e-12  

 alg: bhu
 time update: 0.9177 
 time mult:   0.9177
 error:       0.00000e+00

```

From Matlab one can run the same example using

``>> matlab_example``

```
***************** Bidiag update ******************* 
*         Algorithms: bgu                           
*                     bhu                           
*                                                   
*           Software: Matlab R2023b                     
*             System: MACI64                            
*                                                   
*          Prob size: m=2000                          
*                     n=800                          
*                     nnz(w)=3                     
*                     nnz(p)=2                     
*                                                   
*         Release Aug. 2025                         
*         J.J. Brust (johannesbrust@yahoo.com)      
*************************************************** 

 alg: bgu                         
 time update: 0.0496                
 time mult:   0.0496                
 error:       3.63798e-12                
                                  
 alg: bhu                         
 time update: 0.8650                
 time mult:   0.8650                
 error:       3.63798e-12 

```

**Compilation from source**

*Unix (linux/mac)*

To compile the source files on a linux computer, you can navigate
to folder "include/" and first compile auxiliary subroutines.
The flag "-Ofast" yields aggressive code optimization, and usually the fastest
results but may lead to different floating point behaviour. 
It can be changed to "-O3", or something similar, or even omitted altogether.

    gfortran -c -Ofast givens.f90 kind_parameter.f90

    gfortran -c -Ofast bidiag_up.f90
    gfortran -c -Ofast bidiagh_up.f90

Ensure that a copy of the ".mod" and ".o" files with the sames names
as the ".f90" files of the previous step are in the "test/" directory.
Now you can create an executable test, like for instance

    gfortran -c test_bdup_givapply_qmul_1.f90
    gfortran -o test_bdup test_bdup_givapply_qmul_1.o bidiag_up.o givens.o kind_parameter.o

To run this test you can type

    ./test_bdup

**Compilation of Interfaces**

*Matlab*

64-bit matlab executables for linux and mac are provided as part of this package.
However, you may want (or need) to compile the Matlab interface yourself.
For that the "mex" C compiler has to be setup. The initial step is to 
invoke "mex -setup" in Matlab. Once a C compiler is available (via mex)
in Matlab, you can compile the algorithms like so:

    mex bidiagup.c bidiag_up.o givens.o kind_parameter.o -lgfortran

    mex bhu.c bidiagh_up.o kind_parameter.o -lgfortran

    mex bdup_mulp.c bidiag_up.o givens.o kind_parameter.o -lgfortran

    mex bdup_mulq.c bidiag_up.o givens.o kind_parameter.o -lgfortran

Now you can call the algorithms like

    [BP2,QM,QB,QP2,PM,PP2] = bidiagup(B1,b,c);

    BUP = [diag(B1(n1,n1)),[diag(B1(n1,n1),1);0]];                    
    [B,Y,W] = bhu(BUP,b,c,n);

*Python*

To use the algorithms in Python the `ctypes` library accesses the functions
via a shared object (.so) in linux or a dynamic library (.dylib) on mac.
You will use the object files generated in the section "Compilation from source" from above

Linux

```
gfortran -shared bidiagup.o givens.o kind_parameter.p -o bidiagup.so
gfortran -shared bhu.o kind_parameter.o -o bhu.so
```

MacOS

```
gfortran -shared bidiagup.o givens.o kind_parameter.p -o bidiagup.dylib
gfortran -shared bhu.o kind_parameter.o -o bhu.dylib
```

 ## References
[1] J.J. Brust and M.A. Saunders, Fast and Accurate SVD-Type Updating in Streaming Data, (2025)
