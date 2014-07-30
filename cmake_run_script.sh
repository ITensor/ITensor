#!/bin/bash
#To compile with cmake, make a copy of this file in your build directory and modify it according to your needs.

#Compiler: e.g. icc/icpc, gcc/g++ etc.
export CC=icc
export CXX=icpc

#Set a BLAS/LAPACK vendor: e.g. Generic, Intel, acml, etc...
#(See FindBLAS.cmake for more options.)
#export BLA_VENDOR=Intel
export BLA_VENDOR=Generic

cmake .. -DBLA_VENDOR=$BLA_VENDOR
