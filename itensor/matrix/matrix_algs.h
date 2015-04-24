//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_ALGS__H_
#define __ITENSOR_MATRIX_ALGS__H_

#include "matrix.h"

namespace itensor {

void
diagSymmetric(const matrixref& M,
              matrixref& U,
              vecref& d);

void
diagSymmetric(const matrixref& M,
              matrix& U,
              vec& d);

//orthogonalize the first num columns of m,
//optionally repeating numpass times to reduce roundoff errors
void 
orthog(const matrixref& m, long num, long numpass = 2);

};

#endif
