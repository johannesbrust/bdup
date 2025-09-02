#!/bin/bash
cd ../include
gfortran -c bidiag_up.f90
cp bidiag_up.o ../objs/bidiag_up.o
cp bidiag_up.mod ../test/bidiag_up.mod
cd ../test
gfortran -c test_bdup_givapply_qmul_1.f90
gfortran -o test_bdup_givapply_qmul_1 test_bdup_givapply_qmul_1.o ../objs/bidiag_up.o ../objs/kind_parameter.o ../objs/givens.o
