#include "svd.h"

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
// Making newThresh bigger improves the accuracy but
// makes the algorithm run slower.
//
// If newThresh == 0 the algorithm does only one pass.
//

void
SVD(const MatrixRef& A, Matrix& U, Vector& D, Matrix& V,
    Real newThresh)
    {
    const int n = A.Nrows(), 
              m = A.Ncols();

    if(n > m)
        {
        Matrix At = A.t(), Ut, Vt;
        SVD(At,Vt,D,Ut,newThresh);
        U = Ut.t();
        V = Vt.t();
        return;
        }

    //Form density matrix rho
    Matrix rho = A * A.t();

    //Diagonalize rho
    rho *= -1; //Negative sign sorts evals from > to <
    Vector evals;
    EigenValues(rho,evals,U);

    //Form Vt and fix up its orthogonality
    //(Vt is transpose of V)
    Matrix Vt = A.t() * U;
    Orthog(Vt,n,2); //2 is the number of orthog passes

    //B should be close to diagonal
    //but may not be perfect - fix
    //it up below
    Matrix B = U.t() * A * Vt;

    D = B.Diagonal();
    V = Vt.t();

    if(D(1) == 0 || newThresh == 0) return;

    int start = 2;
    const Real D1 = D(1);
    for(; start < n; ++start)
        {
        if(D(start)/D1 < newThresh) break;
        }

    if(start >= (n-1)) return;

    //
    //Recursively SVD part of B 
    //for greater final accuracy
    //

    Matrix b = B.SubMatrix(start,n,start,n);

    Matrix u,v;
    Vector d;
    SVD(b,u,d,v,newThresh);

    D.SubVector(start,n) = d;

    U.SubMatrix(1,n,start,n) = U.SubMatrix(1,n,start,n) * u;

    V.SubMatrix(start,n,1,m) = v * Vt.t().SubMatrix(start,n,1,m);

    return;
    }
