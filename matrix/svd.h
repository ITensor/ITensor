#ifndef __MATRIXREF_SVD_H
#define __MATRIXREF_SVD_H

#include "matrix.h"

namespace itensor {

class Matrix;
class Vector;

//
// Performs an accurate singular value decomposition
// of a rectangular n x m Matrix A such that
// A = U * D * V
//
// Works by computing U and V such that B = U.t() * A * V.t()
// should be diagonal. In general it won't be after the
// first pass due to errors in EigenValues. So, take the
// part of B that is not diagonal (usually part involving 
// the smallest singular values) and SVD it, etc.
//
// Making newThresh larger improves the accuracy but
// makes the algorithm run slower.
//
// For the special value newThresh == 0 the algorithm does only one pass.
//

void 
SVD(const MatrixRef& A, Matrix& U, Vector& D, Matrix& V,
    Real newThresh = 1E-4);

void 
SVD(const MatrixRef& Are, const MatrixRef& Aim, 
    Matrix& Ure, Matrix& Uim, 
    Vector& D, 
    Matrix& Vre, Matrix& Vim,
    Real newThresh = 1E-4);

void 
SVDComplex(const MatrixRef& Are, const MatrixRef& Aim,
           Matrix& Ure, Matrix& Uim, 
           Vector& D, 
           Matrix& Vre, Matrix& Vim);
           

}

#endif
