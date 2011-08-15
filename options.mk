### User Configurable Options

CCCOM=g++ -m64

PREFIX=$(THIS_DIR)
LIBDIR=$(PREFIX)/lib
INCLUDEDIR=$(PREFIX)/include
BOOST_DIR=$(PREFIX)/../boost

###BLAS/LAPACK Related Options

##For a recent Mac OSX system (include flags intentionally left blank)
PLATFORM=macos
BLAS_LAPACK_LIBFLAGS=-framework Accelerate
BLAS_LAPACK_INCLUDEFLAGS=

##Example using the Intel MKL library
#PLATFORM=intel
#BLAS_LAPACK_INCLUDEFLAGS=-I/sopt/intel/mkl/10.1.0.015/include
#BLAS_LAPACK_LIBFLAGS=-L/sopt/intel/mkl/10.1.0.015/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lgfortran -lpthread
