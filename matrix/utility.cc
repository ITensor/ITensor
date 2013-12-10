// utility.cc -- Various matrix utility routines

#include "matrix.h"
#include "tarray1.h"
#include "minmax.h"
#include <math.h>
#include <fstream>

#include "lapack_wrap.h"

using namespace std;

Real inline
sqr(Real x) { return x*x; }

//
//Orthogonalize num columns of the real matrix M
//
void 
Orthog(const MatrixRef& M,int num,int numpass)
    {
    //const Real tolerance = 1e-10;

    int nkeep;			// Orthogonalize to at most the column dim 
    if (num > M.Nrows() || (num == 0 && M.Ncols() > M.Nrows()))
        _merror("Ncols() > M.Nrows() in Orthog!");

    if (num > 0 && num <= M.Ncols() && num <= M.Nrows())
        nkeep = num;
    else
        nkeep = min(M.Nrows(), M.Ncols());

    Vector dots(nkeep);
    MatrixRef Mcols;
    VectorRef dotsref, coli;
    int i;
    for (i = 1; i <= nkeep; i++)
        {
        coli << M.Column(i);
        Real norm = Norm(coli);
        if(norm == 0.0)
            {
            coli.Randomize();
            norm = Norm(coli);
            }
        coli /= norm;
        if (i == 1)
            continue;

        Mcols << M.Columns(1,i-1);
        dotsref << dots.SubVector(1,i-1);
        int pass;
        for(pass = 1; pass <= numpass; pass++)
            {
            dotsref = Mcols.t() * coli;
            coli -= Mcols * dotsref;
            Real norm = Norm(coli);
            if(norm < 1.0e-3)   // orthog is suspect
            pass--;
            if(norm < 1.0e-10)  // What if a subspace was zero in all vectors?
                {
                coli.Randomize();
                norm = Norm(coli);
                }
            coli /= norm;
            }
        }
    }

//
//Orthogonalize num columns of the complex matrix Mre+i*Mim
//
void 
Orthog(const MatrixRef& Mre, const MatrixRef& Mim,
       int num, int numpass)
    {
    if (num > Mre.Nrows() || (num == 0 && Mre.Ncols() > Mre.Nrows()))
        _merror("Ncols() > M.Nrows() in Orthog!");

#ifdef DEBUG
    if(Mre.Nrows() != Mim.Nrows() || Mre.Ncols() != Mre.Ncols())
        {
        _merror("Mre and Mim must have same dimensions");
        }
#endif

    int nkeep = -1;	// Orthogonalize to at most the column dim 
    if (num > 0 && num <= Mre.Ncols() && num <= Mre.Nrows())
        nkeep = num;
    else
        nkeep = min(Mre.Nrows(), Mim.Ncols());

    Vector redot(nkeep),
           imdot(nkeep);
    MatrixRef Mrecols,
              Mimcols;
    VectorRef dotsre,
              dotsim,
              recol,
              imcol;
    for(int i = 1; i <= nkeep; ++i)
        {
        recol << Mre.Column(i);
        imcol << Mim.Column(i);
        Real norm = sqrt(sqr(Norm(recol))+sqr(Norm(imcol)));
        if(norm == 0.)
            {
            recol.Randomize();
            imcol = 0.;
            norm = Norm(recol);
            }
        recol /= norm;
        imcol /= norm;

        if (i == 1)
            continue;

        Mrecols << Mre.Columns(1,i-1);
        Mimcols << Mim.Columns(1,i-1);
        dotsre << redot.SubVector(1,i-1);
        dotsim << imdot.SubVector(1,i-1);

        for(int pass = 1; pass <= numpass; ++pass)
            {
            dotsre = Mrecols.t() * recol + Mimcols.t() * imcol;
            dotsim = Mrecols.t() * imcol - Mimcols.t() * recol;
            recol -= Mrecols * dotsre - Mimcols * dotsim;
            imcol -= Mrecols * dotsim + Mimcols * dotsre;
            norm = sqrt(sqr(Norm(recol))+sqr(Norm(imcol)));
            if(norm < 1E-3)   // orthog is suspect
                {
                --pass;
                }
            if(norm < 1E-10)  // What if a subspace was zero in all vectors?
                {
                recol.Randomize();
                imcol = 0.;
                norm = Norm(recol);
                }
            recol /= norm;
            imcol /= norm;
            }
        }
    }

void 
QRDecomp(const MatrixRef& M, Matrix& Q, Matrix& R)
    {
    int m = M.Nrows();
    int n = M.Ncols();
    int tlen = min(m,n);
    Vector Tau(tlen); Tau = 0;

    if(m > n)
        {
        Error("Ncols < Nrows, but currently only Ncols >= Nrows case supported");
        }

    Q = M.t();

    int info = 0;

    //Call lapack routine

    dgeqrf_wrapper(&m, &n, Q.Store(), &m, Tau.Store(), &info);

    //int* jpvt = new int[n];
    //for(int i = 0; i < n; ++i) jpvt[i] = 0;
    //dgeqp3_(&m, &n, Q.Store(), &m, jpvt, Tau.Store(), Work.Store(), &lwork, &info);
    //delete jpvt;

    if(info != 0) error("Error in call to dgeqrf_.");

    //Grab R
    R = Matrix(m,n);
    R = 0;
    //Grab elements of R from Q
    for(int i = 1; i <= m; ++i)      
    for(int j = i; j <= n; ++j) 
        {
        R(i,j) = Q(j,i);
        }       

    //Generate Q
    dorgqr_wrapper(&tlen, &tlen, &tlen, Q.Store(), &m, Tau.Store(), &info);
    if(info != 0) error("Error in call to dorgqr_.");

    Q = Q.t().SubMatrix(1,tlen,1,tlen);

    } //void QRDecomp


//
// Eigenvalues and eigenvectors of a real, symmetric matrix A
//
void 
EigenValues(const MatrixRef& A, Vector& D, Matrix& Z)
    {
    LAPACK_INT N = A.Ncols();
    if(N == 0)
      _merror("EigenValues: 0 dimensions matrix");
    if (N != A.Nrows() || A.Nrows() < 1)
	{
	cout << A.Nrows() << " " << A.Ncols() << endl;
	_merror("EigenValues: Input Matrix must be square");
	}

    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT info;
    
    D.ReDimension(N);
    Z = A;

    dsyev_wrapper(&jobz,&uplo,&N,Z.Store(),&N,D.Store(),&info);

    if(info != 0)
        {
        cerr << "info is " << info << endl;
        cout << "info is " << info << endl;
        Error("Got an error code in EigenValues");
        cout << "redoing EigenValues " << endl;
        cerr << "redoing EigenValues " << endl;
        Matrix AA(A);
        for(int i = 1; i <= N; i++)
        for(int j = i+1; j <= N; j++)
            {
            if(AA(i,j) != AA(j,i))
                cout << "Asym: " << i SP j SP AA(i,j) SP AA(j,i) << endl;
            }
        //BackupEigenValues(A,D,Z);
        return;
        }

    //Transpose Z before return
    Z = Z.t();
    }

//
//Compute eigenvalues of arbitrary real matrix A
//
void 
GenEigenValues(const MatrixRef& A, Vector& Re, Vector& Im)
    {
    if(A.Nrows() < 1)
      _merror("GenEigenValues: 0 dimensions matrix");
    if (A.Ncols() != A.Nrows() || A.Nrows() < 1)
        _merror("GenEigenValues: Input Matrix must be square");

    if(A.Ncols() == 1)
        {
        Re = A.Column(1);
        Im = Vector(1); Im(1) = 0.0;
        return;
        }

    LAPACK_INT n = A.Ncols();

    char jobvl = 'N'; //don't compute left evecs
    char jobvr = 'N'; //don't compute right evecs

    Matrix Z = A;
    Re = A.Column(1);
    Im = Re;

    Vector noevecs(2); noevecs = 0; //Shouldn't be referenced by dgeev
    LAPACK_INT info = 0;

    //Call routine
    dgeev_wrapper(&jobvl,&jobvr,&n,Z.Store(),Re.Store(),Im.Store(),noevecs.Store(),noevecs.Store(),&info);
	if(info != 0)
        {
        cerr << "info is " << info << endl;
        _merror("Error condition in dgeev.");
        }
    }

//
//Compute eigenvectors and eigenvalues of arbitrary real matrix A
//
// Convention is that columns of ReV and ImV contain the real and imag
// parts of eigenvectors of A.
//
// If we formally define Q = ReV+i*ImV (where i*i = -1) and 
// D a diagonal matrix containing Re+i*Im on its diagonal, then
// A = Q*D*Q^{-1}
//
// Alternatively A*Q = Q*D
//
void 
GenEigenValues(const MatrixRef& A, Vector& Re, Vector& Im, Matrix& ReV, Matrix& ImV)
    {
    if(A.Nrows() < 1)
      _merror("GenEigenValues: 0 dimensions matrix");
    if (A.Ncols() != A.Nrows() || A.Nrows() < 1)
        _merror("GenEigenValues: Input Matrix must be square");

    /*
    if(A.Ncols() == 1)
        {
        Re = A.Column(1);
        Im = Vector(1); Im(1) = 0.0;
        ReV = A; ReV = 1.;
        ImV = A; ImV = 0.;
        return;
        }
        */

    LAPACK_INT N = A.Ncols();

    char jobvl = 'N'; //don't compute left evecs
    char jobvr = 'V'; //compute right evecs

    Matrix Z = A;
    Z = Z.t();
    Re = A.Column(1);
    Im = Re;

    Vector noLevecs(2); //place holder for left evecs (not computed)
    LAPACK_INT info = 0;

    Matrix evecs(N,N);

    //Call routine
    dgeev_wrapper(&jobvl,&jobvr,&N,Z.Store(),Re.Store(),Im.Store(),noLevecs.Store(),evecs.Store(),&info);
	if(info != 0)
        {
        cout << "info is " << info << endl;
        _merror("Error condition in dgeev.");
        }


    ReV.ReDimension(N,N);
    ImV.ReDimension(N,N);

    evecs = evecs.t();
    //cout << "Evecs = \n" << evecs << endl;

    int n = 1; //nth eigenvalue
    while(n <= N)
        {
        //Check for complex eig pair
        if(Im(n) > 0)
            {
            ReV.Column(n) = evecs.Column(n);
            ReV.Column(n+1) = evecs.Column(n);
            ImV.Column(n) = evecs.Column(n+1);
            ImV.Column(n+1) = -evecs.Column(n+1);
            n += 2;
            }
        else
            {
            ReV.Column(n) = evecs.Column(n);
            ImV.Column(n) = 0.;
            ++n;
            }
        }

    }

//
// Compute eigenvalues and eigenvectors of generalized eigenproblem
// A*x = lambda*B*x
//
void 
GeneralizedEV(const MatrixRef& A, Matrix B, Vector& D, Matrix& Z)
    {
    LAPACK_INT N = A.Ncols();
    if(A.Nrows() < 1)
      _merror("GeneralizedEV: 0 dimensions matrix");
    if (N != A.Nrows() || A.Nrows() < 1)
      _merror("GeneralizedEV: Input Matrix must be square");

    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT info;
    
    D.ReDimension(N);
    Z = A;

    dsygv_wrapper(&jobz,&uplo,&N,Z.Store(),B.Store(),D.Store(),&info);

    if(info != 0)
        {
        cerr << "info is " << info << endl;
        cout << "info is " << info << endl;
        Error("Error in call to dsygv");
        return;
        }

    //Transpose Z before return
    Z = Z.t();
    }

void
ComplexEigenvalues(const MatrixRef& Mre, const MatrixRef& Mim,
                   Vector& revals, Vector& ievals,
                   Matrix& revecs, Matrix& ievecs)
    {
    LAPACK_INT N = Mre.Nrows();
#ifdef DEBUG
    if(Mre.Ncols() != N)
        Error("Mre must be square");
    if(Mim.Nrows() != N || Mim.Ncols() != N)
        Error("Mim must have same dimensions as Mre");
#endif

    Matrix M(N,2*N);
    for(int i = 1; i <= N; ++i)
	for(int j = 1; j <= N; ++j)
        {
	    M(i,2*j-1) = Mre(j,i); 
        M(i,2*j) = Mim(j,i);
        }

    char jobvl = 'N';
    char jobvr = 'V';
    int info = 0;
    
    Vector evals(2*N);
    Matrix evecs(N,2*N);
    LAPACK_COMPLEX l;

    zgeev_wrapper(&jobvl,&jobvr,&N,(LAPACK_COMPLEX*)M.Store(),(LAPACK_COMPLEX*)evals.Store(),&l,
                  (LAPACK_COMPLEX*)evecs.Store(),&info);

    if(info != 0)
        {
        cout << "info is " << info << endl;
        _merror("Error in ComplexEigenvalues");
        }

    revals.ReDimension(N);
    ievals.ReDimension(N);

    revecs.ReDimension(N,N);
    ievecs.ReDimension(N,N);

    for(int i = 1; i <= N; ++i)
        {
        revals(i) = evals(2*i-1);
        ievals(i) = evals(2*i);

        for(int j = 1; j <= N; ++j)
            {
            revecs(j,i) = evecs(i,2*j-1); 
            ievecs(j,i) = evecs(i,2*j);
            }
        }
    }


void 
HermitianEigenvalues(const Matrix& re, const Matrix& im, 
                     Vector& evals,
                     Matrix& revecs, Matrix& ievecs)
    {
    LAPACK_INT N = re.Ncols();
    if(re.Nrows() < 1)
      _merror("HermitianEigenvalues: 0 dimensions re matrix");
    if (N != re.Nrows() || re.Nrows() < 1)
      _merror("HermitianEigenValues: Input Matrix must be square");
    if(im.Ncols() != N || im.Nrows() != N)
      _merror("HermitianEigenValues: im not same dimensions as re");

    Matrix AA(N,2*N);
    for(int i = 1; i <= N; ++i)
	for(int j = 1; j <= N; ++j)
        {
	    AA(i,2*j-1) = re(j,i); 
        AA(i,2*j) = im(j,i);
        }

    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT lwork = max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    LAPACK_COMPLEX work[lwork];
    LAPACK_REAL rwork[lwork];
    LAPACK_INT info;
    
    evals.ReDimension(N);

    zheev_wrapper(&jobz,&uplo,&N,(LAPACK_COMPLEX*)AA.Store(),&N,evals.Store(),work,&lwork,rwork,&info);

    if(info != 0)
        {
        cout << "info is " << info << endl;
        _merror("HermitianEigenvalues: info bad");
        }

    revecs.ReDimension(N,N);
    ievecs.ReDimension(N,N);

    for(int i = 1; i <= N; ++i)
    for(int j = 1; j <= N; ++j)
        {
        revecs(j,i) = AA(i,2*j-1); 
        ievecs(j,i) = AA(i,2*j);
        }

    }

#ifdef I
#undef I
#endif


// Routines needed for Matrix Inverse

// Numerical Recipes (1st. edition) routine for LU decomposition

void ludcmp(Matrix& a,int* indx,Real* d)
    {
    const Real TINY = 1e-20;
    int i, imax=0, j, k;
    Real big, dum, sum, temp;
    Real *vv;
    int n = a.Nrows();

    if (a.Ncols() != n)
	_merror("ludcmp: Must input square matrix");

    Vector V(n);
    vv = V.Store() - 1;
    *d = 1.0;
    for (i = 1; i <= n; i++)
	{
	big = 0.0;
	for (j = 1; j <= n; j++)
	    if ((temp = fabs(a(i, j))) > big)
		big = temp;
	if (big == 0.0)
	    _merror("Singular matrix in routine LUDCMP");
	vv[i] = 1.0 / big;
	}
    for (j = 1; j <= n; j++)
	{
	for (i = 1; i < j; i++)
	    {
// a(i,j) -= a.Row(i).SubVector(1,i-1) * a.Column(j).SubVector(1,i-1);
	    sum = a(i, j);
	    for (k = 1; k < i; k++)
		sum -= a(i, k) * a(k, j);
	    a(i, j) = sum;
	    }
	big = 0.0;
	for (i = j; i <= n; i++)
	    {
	    sum = a(i, j);
	    for (k = 1; k < j; k++)
		sum -= a(i, k) * a(k, j);
	    a(i, j) = sum;
// sum = a.Row(i).SubVector(1,i-1) * a.Column(j).SubVector(1,i-1);
// a(i,j) -= sum;
	    if ((dum = vv[i] * fabs(sum)) >= big)
		{
		big = dum;
		imax = i;
		}
	    }
	if (j != imax)
	    {
	    for (k = 1; k <= n; k++)
		{
		dum = a(imax, k);
		a(imax, k) = a(j, k);
		a(j, k) = dum;
		}
	    *d = -(*d);
	    vv[imax] = vv[j];
	    }
	indx[j] = imax;
	if (a(j, j) == 0.0)
	    a(j, j) = TINY;
	if (j != n)
	    {
	    dum = 1.0 / (a(j, j));
	    for (i = j + 1; i <= n; i++)
		a(i, j) *= dum;
	    }
	}
    }

// NR, 1st edition routine lubksb

void lubksb(Matrix& a,int* indx,Vector& b)
    {
    int i, ii = 0, ip, j;
    Real sum;
    int n = a.Nrows();

    if (a.Ncols() != n)
	_merror("lubksb: must have square input matrix");

    for (i = 1; i <= n; i++)
	{
	ip = indx[i];
	sum = b(ip);
	b(ip) = b(i);
	if (ii)
	    for (j = ii; j <= i - 1; j++)
		sum -= a(i, j) * b(j);
// sum -= a.Row(i).SubVector(ii,i-1) * b.SubVector(ii,i-1);
	else if (sum)
	    ii = i;
	b(i) = sum;
	}
    for (i = n; i >= 1; i--)
	{
	sum = b(i);
	for (j = i + 1; j <= n; j++)
	    sum -= a(i, j) * b(j);
	b(i) = sum / a(i, i);
// b(i) = (b(i) - a.Row(i).SubVector(i+1,n) * b.SubVector(i+1,n)) / a(i,i);
	}
    }

// C++ interfaces to Numerical Recipes routines

Matrix Inverse(const MatrixRef& M)
{
    int n = M.Nrows();
    if(n == 0)
      _merror("Inverse: 0 dimensions matrix");
    if (M.Ncols() != n) 
	{
	cout << "ncols, nrows are " << M.Ncols() SP n << endl;
	_merror("Inverse: Input matrix must be square");
	}

    Matrix a = M;
    Matrix result(n,n);

    int* index = new int[n]; index--;

    Real d;
    Vector col(n);

    ludcmp(a,index,&d);
    
    int j;
    for (j = 1 ; j <= n ; j++ ) {
	col = 0;
	col(j) = 1;
	lubksb(a,index,col);
	result.Column(j) = col;
    }

    delete[] ++index;

    result.MakeTemp();
    return result;
}



void Inverse(const MatrixRef& M, Matrix& result, Real& logdet, Real& sign)
    {
    int n = M.Nrows();
    if (M.Ncols() != n)
	_merror("Inverse: Input matrix must be square");

    Matrix a = M;

    int *index = new int[n]; index--;
    ludcmp(a, index, &sign);

    Real temp; logdet = 0.0;
    int j;
    for (j = 1; j <= n; j++)
	if ((temp = a(j, j)) < 0.0)
	    { sign = -sign; logdet += log(-temp); }
	else
	    logdet += log(temp);

    Vector col(n);
    result.ReDimension(n, n);
    for (j = 1; j <= n; j++)
	{
	col = 0.0;
	col(j) = 1;
	lubksb(a, index, col);
	result.Column(j) = col;
	}
    delete[]++ index;
    }

Matrix Solve(const MatrixRef& M,const MatrixRef& B)
{
    int n = M.Nrows();
    if (M.Ncols() != n || B.Nrows() != n)
      _merror("Solve: Bad input Matrix");

    Matrix a = M;
    Matrix result(n,n);

    int* index = new int[n]; index--;
    Real d;

    ludcmp(a,index,&d);
    Vector col;
    int j;
    for (j = 1 ; j <= n ; j++ ) {
	col = B.Column(j);
	lubksb(a,index,col);
	result.Column(j) = col;
    }

    delete[] ++index;

    result.MakeTemp();
    return result;
}

// Determinant of a Matrix

void Determinant(const MatrixRef& M, Real& logdet, Real& sign)
    {
    int n = M.Nrows();
    if (M.Ncols() != n)
	_merror("Determinant: Input matrix must be square");

    Matrix a = M;
    int *index = new int[n]; index--;

    ludcmp(a, index, &sign);
    Real temp; logdet = 0.0;
    int j;
    for (j = 1; j <= n; j++)
	if ((temp = a(j,j)) < 0.0)
	    { sign = -sign; logdet += log(-temp); }
	else
	    logdet += log(temp);
    delete[] ++index;
    }

Real Determinant(const MatrixRef& M)
    {
    Real logdet,sign;
    Determinant(M,logdet,sign);
    return sign*exp(logdet);    
    }
//#endif      //  -----}

// Set a Matrix of eigenvectors to the unit vectors plus a small random part
// Can be used to get an initial set of vectors for david().

void resetev(MatrixRef& evecs)
    { evecs.Randomize(); evecs *= 0.1; evecs += 1.0; }

// Use heapsort to sort a Vector in place. Adaptation of Num. Rec. 1st edition
// routine sort.c.

void Sort(Vector& V)
    {
    int n = V.Length();
    if (n < 2)
	return;
    Real *ra = V.Store() - 1;
    int l, j, ir, i;
    Real rra;

    l = n / 2 + 1;
    ir = n;
    for (;;)
	{
	if (l > 1)
	    rra = ra[--l];
	else
	    {
	    rra = ra[ir];
	    ra[ir] = ra[1];
	    if (--ir == 1)
		{
		ra[1] = rra;
		return;
		}
	    }
	i = l;
	j = l + l;
	while (j <= ir)
	    {
	    if (j < ir && ra[j] < ra[j + 1])
		++j;
	    if (rra < ra[j])
		{
		ra[i] = ra[j];
		j += (i = j);
		}
	    else
		j = ir + 1;
	    }
	ra[i] = rra;
	}
    }

// Sort two arrays, a Vector and an IntArray1

void Sort(Vector& V, IntArray1& ind)
    {
    int n = V.Length();
    Real *ra = V.Store() - 1;
    int *rb = (int *) ind.Store() - 1;
    int l, j, ir, i;
    Real rra;
    int rrb;

    if (n != ind.Size())
	_merror("Sort(Vector,IntArray1: bad input sizes");

    if (n < 2)
	return;
    l = (n >> 1) + 1;
    ir = n;
    for (;;)
	{
	if (l > 1)
	    {
	    rra = ra[--l];
	    rrb = rb[l];
	    }
	else
	    {
	    rra = ra[ir];
	    rrb = rb[ir];
	    ra[ir] = ra[1];
	    rb[ir] = rb[1];
	    if (--ir == 1)
		{
		ra[1] = rra;
		rb[1] = rrb;
		return;
		}
	    }
	i = l;
	j = l << 1;
	while (j <= ir)
	    {
	    if (j < ir && ra[j] < ra[j + 1])
		++j;
	    if (rra < ra[j])
		{
		ra[i] = ra[j];
		rb[i] = rb[j];
		j += (i = j);
		}
	    else
		j = ir + 1;
	    }
	ra[i] = rra;
	rb[i] = rrb;
	}
    }

/*

int FFT(const VectorRef& in, Vector& outre, Vector& outim)
    {
    static Vector intemp,outtemp;
    intemp.ReduceDimension(in.Length());
    outtemp.ReduceDimension(in.Length());
    intemp = in;
    int status, n = in.Length(), n2=n/2, stride = 1;
    char ifmt = 'R', ofmt = 'R', direction = 'F';
    status = dfft_(&ifmt,&ofmt,&direction,intemp.Store(),
			outtemp.Store(),&n,&stride);
    outre = outtemp.SubVector(1,n2+1);
    outim.ReduceDimension(n2+1);
    outim[0] = 0.0; outim[n2] = 0.0;
    int i;
    for(i=1; i < n2; i++)
	outim[i] = outtemp[n-i];
    return status;
    }

#endif
*/

Matrix Exp(const MatrixRef& M)
    {
    int n = M.Nrows();
    Matrix evec(n, n);
    Matrix y = M;
    Vector eval(n);
    EigenValues(y, eval, evec);
    int i;
    for (i = 1; i <= n; i++)
        y.Column(i) = evec.Column(i) * exp(eval(i));
    return y * Transpose(evec);
    }



