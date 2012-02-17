#ifndef _SVD_H
#define _SVD_H

#include "matrix.h"

//
// Performs an accurate singular value decomposition
// of a rectangular n x m Matrix A such that
// A = U * D * V
//
// Making newThresh bigger improves the accuracy but
// makes the algorithm run slower.
//
// If newThresh == 0 the algorithm does only one pass.
//

void
SVD(const Matrix& A, Matrix& U, Vector& D, Matrix& V,
    Real newThresh = 1E-4);

#endif
