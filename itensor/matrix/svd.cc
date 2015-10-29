#include "svd.h"

#include <fstream>

#include "itensor/matrix/lapack_wrap.h"
//#include "itensor/util/print.h"

using namespace std;

namespace itensor {

void
writeMatrix(const string& fname, const Matrix& M)
    {
    std::ofstream s(fname.c_str());

    int nr = M.Nrows();
    int nc = M.Ncols();
    s.write((char*)&nr,sizeof(nr));
    s.write((char*)&nc,sizeof(nc));

    Real val;
    for(int j = 1; j <= nr; ++j)
    for(int k = 1; k <= nc; ++k)
        {
        val = M(j,k);
        s.write((char*)&val,sizeof(val));
        }
    s.close();
    }

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
checksvd(const MatrixRef& A, Matrix& U, Vector& D, Matrix& V)
    {
    Matrix Ach = U;
    for(int i = 1; i <= D.Length(); i++)
	Ach.Column(i) *= D(i);
    Ach = Ach * V;
    Ach -= A;
    Real nor = Norm(A.TreatAsVector());
    if(nor != 0.0)
        cout << "relative error with sqrt in low level svd is " << Norm(Ach.TreatAsVector())/nor << endl;
    }

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
#ifdef CHKSVD
        checksvd(A,U,D,V);
#endif
        return;
        }

    //Form 'density matrix' rho
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

    if(D(1) == 0 || newThresh == 0)
        {
#ifdef CHKSVD
        checksvd(A,U,D,V);
#endif
        return;
        }

    int start = 2;
    const Real D1 = D(1);
    for(; start < n; ++start)
        {
        if(D(start)/D1 < newThresh) break;
        }

    if(start >= (n-1)) 
        {
#ifdef CHKSVD
        checksvd(A,U,D,V);
#endif
        return;
        }

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

    V.SubMatrix(start,n,1,m) = v * V.SubMatrix(start,n,1,m);

#ifdef CHKSVD
	checksvd(A,U,D,V);
#endif

    return;
    }

void 
SVD(const MatrixRef& Are, const MatrixRef& Aim, 
    Matrix& Ure, Matrix& Uim, 
    Vector& D, 
    Matrix& Vre, Matrix& Vim,
    Real newThresh)
    {
    const int n = Are.Nrows(), 
              m = Are.Ncols();

#ifdef DEBUG
    if(Aim.Nrows() != n || Aim.Ncols() != m)
        {
        Error("Aim must have same dimensions as Are");
        }
#endif

    if(n > m)
        {
        Matrix Aret = Are.t(),
               Aimt = -Aim.t(),
               Uret, Uimt,
               Vret, Vimt;
        SVD(Aret,Aimt,Vret,Vimt,D,Uret,Uimt,newThresh);
        Ure = Uret.t();
        Uim = -Uimt.t();
        Vre = Vret.t();
        Vim = -Vimt.t();
        return;
        }

    //Form 'density matrix' rho
    Matrix rhore = Are * Are.t() + Aim * Aim.t(),
           rhoim = Aim * Are.t() - Are * Aim.t();


    //Diagonalize rho
    rhore *= -1; //Negative sign sorts evals from > to <
    rhoim *= -1;
    Vector evals;
    HermitianEigenvalues(rhore,rhoim,evals,Ure,Uim);

    //Form Vt and fix up its orthogonality
    //(Vt is transpose of V)
    Matrix Vret = Are.t() * Ure + Aim.t() * Uim,
           Vimt = -Are.t() * Uim + Aim.t() * Ure;
    Orthog(Vret,Vimt,n,2); //2 is the number of orthog passes

    //B should be close to diagonal
    //but may not be perfect - fix it up below
    Matrix Bre = Ure.t()*Are*Vret + Ure.t()*Aim*Vimt + Uim.t()*Aim*Vret - Uim.t()*Are*Vimt,
           Bim = Ure.t()*Aim*Vret - Ure.t()*Are*Vimt - Uim.t()*Are*Vret - Uim.t()*Aim*Vimt;

    Vre = Vret.t();
    Vim = Vimt.t();

    D = Bre.Diagonal();

    if(D(1) == 0 || newThresh == 0)
        {
        return;
        }

    int start = 2;
    const Real D1 = D(1);
    for(; start < n; ++start)
        {
        if(D(start)/D1 < newThresh) break;
        }

    if(start >= (n-1)) 
        {
        return;
        }

    //
    //Recursively SVD part of B 
    //for greater final accuracy
    //

    Matrix bre = Bre.SubMatrix(start,n,start,n),
           bim = Bim.SubMatrix(start,n,start,n);

    Matrix ure,uim,
           vre,vim;
    Vector d;
    SVD(bre,bim,ure,uim,d,vre,vim,newThresh);

    D.SubVector(start,n) = d;

    //Need to copy U and V real part since first of each
    //pair of lines below modified real part, which
    //is used again
    Matrix pUre(Ure),
           pVre(Vre);

    Ure.SubMatrix(1,n,start,n) = pUre.SubMatrix(1,n,start,n) * ure - Uim.SubMatrix(1,n,start,n) * uim;
    Uim.SubMatrix(1,n,start,n) = Uim.SubMatrix(1,n,start,n) * ure + pUre.SubMatrix(1,n,start,n) * uim;

    Vre.SubMatrix(start,n,1,m) = vre * pVre.SubMatrix(start,n,1,m) - vim * Vim.SubMatrix(start,n,1,m);
    Vim.SubMatrix(start,n,1,m) = vim * pVre.SubMatrix(start,n,1,m) + vre * Vim.SubMatrix(start,n,1,m);

    return;
    }


void
SVDComplex(const MatrixRef& Are, const MatrixRef& Aim, 
           Matrix& Ure, Matrix& Uim, 
           Vector& d, 
           Matrix& Vre, Matrix& Vim)
    {
    LAPACK_INT m = Are.Nrows(), 
               n = Are.Ncols(); 
#ifdef DEBUG
    if(Aim.Nrows() != m || Aim.Ncols() != n)
        {
        Error("Aim must have same dimensions as Are");
        }
#endif

    if(m < n)
        {
        Matrix mret = Are.t(), 
               mimt = -Aim.t(),
               UUre,UUim,
               VVre,VVim;
        SVDComplex(mret,mimt, UUre, UUim, d, VVre, VVim);
        Vre = UUre.t();
        Vim = -UUim.t();
        Ure = VVre.t();
        Uim = -VVim.t();
        return;
        }

    char jobz = 'S';

    Matrix AA(n,2*m);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
        {
	    AA(i,2*j-1) = Are(j,i); 
        AA(i,2*j) = Aim(j,i);
        }

    Matrix UU(n,2*m), VV(n,2*n);
    d.ReDimension(n);
    LAPACK_INT info = 0;

    zgesdd_wrapper(&jobz,&m,&n,
                   (LAPACK_COMPLEX*)AA.Store(),
                   d.Store(), 
                   (LAPACK_COMPLEX*)UU.Store(),
                   (LAPACK_COMPLEX*)VV.Store(),
                   &info);

    if(info != 0) 
        {
        cout << "info = " << info << endl;
        Error("Error condition in zgesdd");
        }

    Ure.ReDimension(m,n);
    Uim.ReDimension(m,n);
    Vre.ReDimension(n,n);
    Vim.ReDimension(n,n);

    for(int i = 1; i <= n; ++i)
    for(int j = 1; j <= m; ++j)
        {
        Ure(j,i) = UU(i,2*j-1); 
        Uim(j,i) = UU(i,2*j);
        }

    for(int i = 1; i <= n; ++i)
    for(int j = 1; j <= n; ++j)
        {
        Vre(j,i) = VV(i,2*j-1); 
        Vim(j,i) = VV(i,2*j);
        }
    }

}
