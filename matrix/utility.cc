// utility.cc -- Various matrix utility routines

#include "matrix.h"
#include "tarray1.h"
#include "minmax.h"
#include <math.h>
#include <fstream>

#include "lapack_wrap.h"

using namespace std;

void Orthog(const MatrixRef& M,int num,int numpass)	// Orthonormalize a Matrix M to num cols
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

void 
QRDecomp(const MatrixRef& M, Matrix& Q, Matrix& R)
    {
    int m = M.Nrows();
    int n = M.Ncols();
    int tlen = min(m,n);
    Vector Tau(tlen); Tau = 0;

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
    R = Matrix(tlen,tlen);
    R = 0;
    //Grab elements of R from Q
    for(int i = 1; i <= tlen; ++i)      
    for(int j = i; j <= tlen; ++j) 
        {
        R(i,j) = Q(j,i);
        }       

    //Generate Q
    dorgqr_wrapper(&m, &n, &tlen, Q.Store(), &m, Tau.Store(), &info);
    if(info != 0) error("Error in call to dorgqr_.");

    Q = Q.t();

    } //void QRDecomp

void BackupEigenValues(const MatrixRef& A, Vector& D, Matrix& Z);

void 
EigenValues(const MatrixRef& A, Vector& D, Matrix& Z)
    {
    LAPACK_INT N = A.Ncols();
    if (N != A.Nrows() || A.Nrows() < 1)
      _merror("EigenValues: Input Matrix must be square");

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
        cout << "redoing EigenValues " << endl;
        cerr << "redoing EigenValues " << endl;
        Matrix AA(A);
        for(int i = 1; i <= N; i++)
        for(int j = i+1; j <= N; j++)
            {
            if(AA(i,j) != AA(j,i))
                cout << "Asym: " << i SP j SP AA(i,j) SP AA(j,i) << endl;
            }
        BackupEigenValues(A,D,Z);
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
// Compute eigenvalues and eigenvectors of generalized eigenproblem
// A*x = lambda*B*x
//
void 
GeneralizedEV(const MatrixRef& A, Matrix B, Vector& D, Matrix& Z)
    {
    LAPACK_INT N = A.Ncols();
    if (N != A.Nrows() || A.Nrows() < 1)
      _merror("EigenValues: Input Matrix must be square");

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


#include <complex>
#include <vector>
typedef std::complex<double> Complex;
static Complex I(0.0,1.0), C1(1.0,0.0),C0(0.0,0.0);
inline Real sqr(Real a) { return a*a; }

class ComplexVector
    {
public:
    std::vector<Complex> dat;
    ComplexVector(int n=1) : dat(n,C0) { }
    Complex& operator()(int i) { return dat[i-1]; }			// access 1 ... n
    Complex operator()(int i) const { return dat[i-1]; }
    int Length() const { return dat.size(); }
    Vector RealVec()
	{
	int n = Length(); Vector re(n);
	for(int i = 0; i < n; i++)
	    re.el(i) = real(dat[i]);
	return re;
	}
    Vector ImVec()
	{
	int n = Length(); Vector im(n);
	for(int i = 0; i < n; i++)
	    im.el(i) = imag(dat[i]);
	return im;
	}
    Complex operator*(const ComplexVector& other)		// applies conj to first vec
	{
	Complex res = 0;
	int n = Length();
	for(int i = 0; i < n; i++)
	    res += conj(dat[i]) * other.dat[i];
	return res;
	}
    };

class ComplexMatrix;
class CMHelper
    {
public:
    ComplexMatrix *p;
    int r;
    CMHelper(ComplexMatrix *pp, int rr) : p(pp), r(rr) {}
    inline Complex& operator[](int c);
    };

class ComplexMatrix
    {
public:
    std::vector<Complex> dat;
    int nrow, ncol;
    int index(int r, int c) const { 
	if(r > nrow || c > ncol) error("bac index");
	return c + (r-1)*ncol; }

    ComplexMatrix(int nr=1, int nc=1) : dat(nr*nc,C0), nrow(nr), ncol(nc) { }

    ComplexMatrix(const Matrix& re, const Matrix& im) : dat(re.Nrows()*re.Ncols())
	{
	nrow = re.Nrows();
	ncol = re.Ncols();
	for(int i = 1; i <= nrow; i++)
	    for(int j = 1; j <= ncol; j++)
		(*this)(i,j) = Complex(re(i,j),im(i,j));
	}

    Complex& operator()(int r, int c) { return dat[index(r,c)-1]; }
    Complex operator()(int r, int c) const { return dat[index(r,c)-1]; }
    CMHelper operator[](int r) { return CMHelper(this,r); }
    Complex& el(int r, int c) { return dat[index(r+1,c+1)-1]; }
    Complex el(int r, int c) const { return dat[index(r+1,c+1)-1]; }

    Matrix RealMat()
	{
	Matrix re(nrow,ncol);
	for(int i = 1; i <= nrow; i++)
	    for(int j = 1; j <= ncol; j++)
		re(i,j) = real(dat[index(i,j)-1]);
	return re;
	}
    Matrix ImMat()
	{
	Matrix im(nrow,ncol);
	for(int i = 1; i <= nrow; i++)
	    for(int j = 1; j <= ncol; j++)
		im(i,j) = imag(dat[index(i,j)-1]);
	return im;
	}
    ComplexVector operator*(const ComplexVector& v)
	{
	ComplexMatrix &This(*this);
	ComplexVector res(nrow);
	for(int i = 1; i <= nrow; i++)
	    for(int j = 1; j <= ncol; j++)
		res(i) += This(i,j) * v(j);
	return res;
	}
    ComplexMatrix Inverse()
	{
	if(nrow != ncol) error("bad nrow ncol in Inverse");
	Matrix re = RealMat(), im = ImMat();
	Matrix big(2*nrow,2*nrow);
	big.SubMatrix(1,nrow,1,nrow) = re;
	big.SubMatrix(nrow+1,2*nrow,nrow+1,2*nrow) = re;
	big.SubMatrix(1,nrow,nrow+1,2*nrow) = im;
	big.SubMatrix(nrow+1,2*nrow,1,nrow) = -im;
	Matrix ibig = ::Inverse(big);
	re = ibig.SubMatrix(1,nrow,1,nrow);
	im = ibig.SubMatrix(1,nrow,nrow+1,2*nrow);
	ComplexMatrix res(nrow,nrow);
	for(int i = 1; i <= nrow; i++)
	    for(int j = 1; j <= ncol; j++)
		res(i,j) = Complex(re(i,j),im(i,j));
	return res;
	}
    };

inline Complex& CMHelper::operator[](int c) { return (*p)(r+1,c+1); }

void 
HermitianEigenvalues(const Matrix& re, const Matrix& im, Vector& evals,
                     Matrix& revecs, Matrix& ievecs)
    {
    LAPACK_INT N = re.Ncols();
    if (N != re.Nrows() || re.Nrows() < 1)
      _merror("HermitianEigenValues: Input Matrix must be square");
    if(im.Ncols() != N || im.Nrows() != N)
      _merror("HermitianEigenValues: im not same dimensions as re");

    Matrix imt(im.t());
    ComplexMatrix H(re,imt), evecs(re,imt);

    char jobz = 'V';
    char uplo = 'U';
    LAPACK_INT lwork = max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    LAPACK_COMPLEX work[lwork];
    LAPACK_REAL rwork[lwork];
    LAPACK_INT info;
    
    evals.ReDimension(N);

    zheev_wrapper(&jobz,&uplo,&N,(LAPACK_COMPLEX*)&(H.dat[0]),&N,evals.Store(),work,&lwork,rwork,&info);
    revecs = H.RealMat().t();
    ievecs = H.ImMat().t();

    if(info != 0)
        {
        cerr << "info is " << info << endl;
        cout << "info is " << info << endl;
        //cout << "redoing EigenValues " << endl;
        //cerr << "redoing EigenValues " << endl;
        //Matrix AA(A);
        //for(int i = 1; i <= N; i++)
        //for(int j = i+1; j <= N; j++)
        //if(AA(i,j) != AA(j,i))
            //cout << "Asym: " << i SP j SP AA(i,j) SP AA(j,i) << endl;
        //BackupEigenValues(A,D,Z);
        //return;
        _merror("EigenValues: info bad");
        }
    }

// Possibly extra stuff for EigenValues:
// This is grabbed from evalue.cc:

#include "precisio.h"

inline Real square(Real x) { return x*x; }
inline Real sign(Real x, Real y)
   { return (y>=0) ? x : -x; }                    // assume x >=0

typedef unsigned char Boolean;
const Boolean True = 1;
const Boolean False = 0;

double dotprod(double *a,double *b, int l);
double dotprod2(double *a,int inca, double *b,int incb,  int l);
    
static void tred2(const MatrixRef& A, Vector& D, Vector& E, Matrix& Z)
    {
    Real tol = 1e-200;
//    FloatingPointPrecision::Minimum () / FloatingPointPrecision::Epsilon ();
    int n = A.Nrows ();		// Z.ReDimension(n,n); Z.Inject(A); 
    Z = A;
    D.ReDimension (n);
    E.ReDimension (n);
    Real *z = Z.Store ();

    int i;
    for (i = n - 1; i > 0; i--)	// i=0 is excluded 
	{
	Real f = Z.el (i, i - 1);
	Real g = 0.0;
	int k = i - 1;
	Real *zik = z + i * n;
	while (k--)
	    g += square (*zik++);
	Real h = g + square (f);
	if (g <= tol)
	    {
	    E.el (i) = f;
	    h = 0.0;
	    }
	else
	    {
	    g = sign (-sqrt (h), f);
	    E.el (i) = g;
	    h -= f * g;
	    Z.el (i, i - 1) = f - g;
	    f = 0.0;
	    Real *zji = z + i;
	    Real *zij = z + i * n;
	    Real *ej = E.Store ();
	    int j;
	    for (j = 0; j < i; j++)
		{
		*zji = (*zij++) / h;
		Real *zjk = z + j * n;
		zik = z + i * n;
		g = dotprod(zjk,zik,j);
		zjk += j;
		zik += j;
		g += dotprod2(zjk,n,zik,1,i-j);
		*ej++ = g / h;
		f += g * (*zji);
		zji += n;
		}
	    Real hh = f / (h + h);
	    zij = z + i * n;
	    ej = E.Store ();
	    for (j = 0; j < i; j++)
		{
		f = *zij++;
		g = *ej - hh * f;
		*ej++ = g;
		Real *zjk = z + j * n;
		Real *zik = z + i * n;
		Real *ek = E.Store ();
		k = j + 1;
		int kk = k & 1;
		k -= kk;
		while (kk--)
		    {
		    Real etemp = *ek;
		    Real ziktemp = *zik;
		    Real zjktemp = *zjk;
		    Real fek = f * etemp; 
		    Real gzik = g * ziktemp; ek++; zik++;
		    zjktemp -= fek;
		    zjktemp -= gzik;
		    *zjk = zjktemp;
		    zjk++;
		    // *zjk++ -= (f * (*ek++) + g * (*zik++));
		    }
		while (k > 0)
		    {
		    k -= 2;
		    Real etemp = *ek;
		    Real etemp2 = *(ek+1);
		    Real ziktemp = *zik; 
		    Real ziktemp2 = *(zik+1);
		    Real zjktemp = *zjk;
		    Real zjktemp2 = *(zjk+1);
		    Real fek = f * etemp; ek += 2; zik += 2;
		    Real fek2 = f * etemp2; zjk += 2;
		    Real gzik = g * ziktemp;
		    Real gzik2 = g * ziktemp2;
		    zjktemp -= fek;
		    zjktemp2 -= fek2;
		    zjktemp -= gzik;
		    zjktemp2 -= gzik2;
		    *(zjk-2) = zjktemp;
		    *(zjk-1) = zjktemp2;
		    
		    // *zjk++ -= (f * (*ek++) + g * (*zik++));
		    }
		}
	    }
	D.el (i) = h;
	}

    D.el (0) = 0.0;
    E.el (0) = 0.0;
    for (i = 0; i < n; i++)
	{
	if (D.el (i) != 0.0)
	    {
	    for (int j = 0; j < i; j++)
		{
		Real *zik = z + i * n;
		Real *zkj = z + j;
		Real g = dotprod2(zik,1,zkj,n,i);
		void daxpy(int n,double a,double *x, int incx, 
					double *y,int incy);
		daxpy(i,-g,z+i,n,z+j,n);
		}
	    }
	Real *zij = z + i * n;
	Real *zji = z + i;
	int j = i;
	while (j--)
	    {
	    *zij++ = 0.0;
	    *zji = 0.0;
	    zji += n;
	    }
	D.el (i) = *zij;
	*zij = 1.0;
	}
    }

static void dotranspose(Real *z, int n)
    {
    int i,j;
    for(i = 0; i < n; i++)
	{
	Real *zinj = z+i*n+i+1;
	Real *zjni = z+(i+1)*n+i;
	for(j = i+1; j < n; j++)
	    {
	    // Real p = z[i*n+j];
	    // z[i*n+j] = z[j*n+i];
	    // z[j*n+i] = p;
	    Real p = *zinj;
	    Real q = *zjni;
	    *zjni = p;
	    *zinj = q;
	    zinj++; zjni += n;
	    }
	}
    }
static void tql2(Vector& D, Vector& E, Matrix& Z)
    {
//   Tracer et("Evalue(tql2)"); 
    Real eps = FloatingPointPrecision::Epsilon ();
    int n = D.Length ();
    Real *z = Z.Store ();
    dotranspose(z,n);
    int l;
    for (l = 1; l < n; l++)
	E.el (l - 1) = E.el (l);
    Real b = 0.0;
    Real f = 0.0;
    E.el (n - 1) = 0.0;
    for (l = 0; l < n; l++)
	{
	int i, j;
	Real & dl = D.el (l);
	Real & el = E.el (l);
	Real h = eps * (fabs (dl) + fabs (el));
	if (b < h)
	    b = h;
	int m;
	for (m = l; m < n; m++)
	    if (fabs (E.el (m)) <= b)
		break;
	Boolean test = False;
	for (j = 0; j < 30; j++)
	    {
	    if (m == l)
		{
		test = True;
		break;
		}
	    Real & dl1 = D.el (l + 1);
	    Real g = dl;
	    Real p = (dl1 - g) / (2.0 * el);
	    Real r = sqrt (p * p + 1.0);
	    dl = el / (p < 0.0 ? p - r : p + r);
	    Real h = g - dl;
	    f += h;
	    Real *dlx = &dl1;
	    i = n - l - 1;
	    while (i--)
		*dlx++ -= h;

	    p = D.el (m);
	    Real c = 1.0;
	    Real s = 0.0;
	    for (i = m - 1; i >= l; i--)
		{
		Real ei = E.el (i);
		Real di = D.el (i);
		Real & ei1 = E.el (i + 1);
		g = c * ei;
		h = c * p;
		if (fabs (p) >= fabs (ei))
		    {
		    c = ei / p;
		    r = sqrt (c * c + 1.0);
		    ei1 = s * p * r;
		    s = c / r;
		    c = 1.0 / r;
		    }
		else
		    {
		    c = p / ei;
		    r = sqrt (c * c + 1.0);
		    ei1 = s * ei * r;
		    s = 1.0 / r;
		    c /= r;
		    }
		p = c * di - s * g;
		D.el (i + 1) = h + s * (c * g + s * di);

		/*
		Real *zki = z + i;
		Real *zki1 = zki + 1;
		int k = n;
		while (k--)
		    {
		    h = *zki1;
		    *zki1 = s * (*zki) + c * h;
		    *zki = c * (*zki) - s * h;
		    zki += n;
		    zki1 += n;
		    }
		*/
		void rotate22(Real *, Real *,Real,Real,int);
		Real *zki = z + i*n;
		Real *zki1 = zki + n;
		rotate22(zki,zki1,c,s,n);
		/*
		int k = n;
		while (k--)
		    {
		    h = *zki1;
		    *zki1 = s * (*zki) + c * h;
		    *zki = c * (*zki) - s * h;
		    zki++;
		    zki1++;
		    }
		*/
		}
	    el = s * p;
	    dl = c * p;
	    if (fabs (el) <= b)
		{
		test = True;
		break;
		}
	    }
	if (!test)
	    _merror ("tql2:Convergence Error");
	dl += f;
	}

    int i;
    for (i = 0; i < n; i++)
	{
	int k = i;
	Real p = D.el (i);
	int j;
	for (j = i + 1; j < n; j++)
	    {
	    if (D.el (j) < p)
		{
		k = j;
		p = D.el (j);
		}
	    }
	if (k != i)
	    {
	    D.el (k) = D.el (i);
	    D.el (i) = p;
	    int j = n;
	    /*
	    Real *zji = z + i;
	    Real *zjk = z + k;
	    while (j--)
		{
		p = *zji;
		*zji = *zjk;
		*zjk = p;
		zji += n;
		zjk += n;
		}
	    */
	    Real *zji = z + i*n;
	    Real *zjk = z + k*n;
	    while (j--)
		{
		p = *zji;
		*zji = *zjk;
		*zjk = p;
		zji++;
		zjk++;
		}
	    }
	}

    dotranspose(z,n);
    }

/*
static void tred3(const Matrix& X, Vector& D, Vector& E, Matrix& A)
{
   Real tol =
      FloatingPointPrecision::Minimum()/FloatingPointPrecision::Epsilon();
   int n = X.Nrows(); A = X; D.ReDimension(n); E.ReDimension(n);
   Real* ei = E.Store() + n;
   int i;
   for (i = n-1; i >= 0; i--)
   {
      Real h = 0.0; Real f;
      Real* d = D.Store(); Real* a = A.Store() + (i*(i+1))/2; int k = i;
      while (k--) { f = *a++; *d++ = f; h += square(f); }
      if (h <= tol) { *(--ei) = 0.0; h = 0.0; }
      else
      {
	 Real g = sign(-sqrt(h), f); *(--ei) = g; h -= f*g;
         f -= g; *(d-1) = f; *(a-1) = f; f = 0.0;
         Real* dj = D.Store(); Real* ej = E.Store();
	 int j;
         for (j = 0; j < i; j++)
         {
            Real* dk = D.Store(); Real* ak = A.Store()+(j*(j+1))/2;
            Real g = 0.0; k = j;
            while (k--)  g += *ak++ * *dk++;
            k = i-j; int l = j; 
            while (k--) { g += *ak * *dk++; ak += ++l; }
	    g /= h; *ej++ = g; f += g * *dj++;
         }  
	 Real hh = f / (2 * h); Real* ak = A.Store();
         dj = D.Store(); ej = E.Store();
         for (j = 0; j < i; j++)
         {
	    f = *dj++; g = *ej - hh * f; *ej++ = g;
            Real* dk = D.Store(); Real* ek = E.Store(); k = j+1;
	    while (k--) { *ak++ -= (f * *ek++ + g * *dk++); }
	 }
      }
      *d = *a; *a = h;
   }
}
*/

/*
static void tql1(Vector& D, Vector& E)
{
//   Tracer et("Evalue(tql1)");
   Real eps = FloatingPointPrecision::Epsilon();
   int n = D.Length();
   int l;
   for (l=1; l<n; l++) E.el(l-1) = E.el(l);
   Real b = 0.0; Real f = 0.0; E.el(n-1) = 0.0;
   for (l=0; l<n; l++)
   {
      int i,j;
      Real& dl = D.el(l); Real& el = E.el(l);
      Real h = eps * ( fabs(dl) + fabs(el) );
      if (b < h) b = h;
      int m;
      for (m=l; m<n; m++) if (fabs(E.el(m)) <= b) break;
      Boolean test = False;
      for (j=0; j<30; j++)
      {
         if (m==l) { test = True; break; }
         Real& dl1 = D.el(l+1);
	 Real g = dl; Real p = (dl1-g) / (2.0*el); Real r = sqrt(p*p + 1.0);
	 dl = el / (p < 0.0 ? p-r : p+r); Real h = g - dl; f += h;
         Real* dlx = &dl1; i = n-l-1; while (i--) *dlx++ -= h;

	 p = D.el(m); Real c = 1.0; Real s = 0.0;
	 for (i=m-1; i>=l; i--)
	 {
            Real ei = E.el(i); Real di = D.el(i);
            Real& ei1 = E.el(i+1);
	    g = c * ei; h = c * p;
	    if ( fabs(p) >= fabs(ei))
	    {
	       c = ei / p; r = sqrt(c*c + 1.0); 
               ei1 = s*p*r; s = c/r; c = 1.0/r;
	    }
	    else
	    {
	       c = p / ei; r = sqrt(c*c + 1.0);
	       ei1 = s * ei * r; s = 1.0/r; c /= r;
	    }
	    p = c * di - s*g; D.el(i+1) = h + s * (c*g + s*di);
	 }
	 el = s*p; dl = c*p;
	 if (fabs(el) <= b) { test = True; break; }
      }
      if (!test) _merror("tql1:Convergence error");
      Real p = dl + f;
      test = False;
      for (i=l; i>0; i--)
      {
         if (p < D.el(i-1)) D.el(i) = D.el(i-1);
         else { test = True; break; }
      }
      if (!test) i=0;
      D.el(i) = p;
   }
}
*/

void BackupEigenValues(const MatrixRef& A, Vector& D, Matrix& Z)
{
    if (A.Ncols() != A.Nrows() || A.Nrows() < 1)
      _merror("BackupEigenValues: Input Matrix must be square");
    Vector E; tred2(A, D, E, Z); tql2(D, E, Z);
}

// Routines needed for Matrix Inverse

// Numerical Recipes (1st. edition) routine for LU decomposition

void ludcmp(Matrix& a,int* indx,Real* d)
    {
    const Real TINY = 1e-20;
    int i, imax, j, k;
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
    if (M.Ncols() != n) _merror("Inverse: Input matrix must be square");

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





void rotate22(double *zki,double *zki1,double c,double s,int n)
    {
    register double h1,h2,h3,h4,t1,t2,t3,t4,r1,r2,r3,r4,s1,s2,s3,s4;
    register double q1,q2,q3,q4,z1,z2,z3,z4;
    double *pzki1 = zki1 - 1;
    double *pzki = zki - 1;
    static double junk[10];
    int nr = n&3 + 4;	// Same as n%4 + 4
    int nmain = n-nr;
    if(nmain < 0) nmain = 0;
    int k;
    if(nmain > 0)
	{
	h1 = *zki1; t1 = *zki;
	zki1++; zki++;
	h2 = *zki1; t2 = *zki;
	zki1++; zki++;
	h3 = *zki1; t3 = *zki;
	zki1++; zki++;
	h4 = *zki1; 
	zki1++;
	k = nmain;
	s1 = s4 = z4 = z2 = z3 = z1 = 0.0;
	pzki = junk;
	pzki1 = junk+4;
	do  {
	    r1 = c*h1;			t4 = *zki;	k -= 4;
	    r2 = s*h1;			zki++;       	
	    r3 = s*t1;	z1 += z3;	h1 = *zki1;     *(pzki1+1) = s1;	
	    r4 = c*t1;	z4 -= z2;	zki1++;      	*(pzki+1) = s4;     
	    q1 = c*h2;			t1 = *zki;      
	    q2 = s*h2;			zki++;       	
	    q3 = s*t2;	r1 += r3;	h2 = *zki1;	*(pzki1+2) = z1;	
	    q4 = c*t2;	r4 -= r2;	zki1++;      	*(pzki+2) = z4;  	
	    s1 = c*h3;			t2 = *zki;	pzki1 = zki1 - 6;
	    s2 = s*h3;			zki++; 		pzki = zki - 6;
	    s3 = s*t3;	q1 += q3;	h3 = *zki1;	*pzki1 = r1;	
	    s4 = c*t3;	q4 -= q2;	zki1+= 2;	*pzki = r4;
	    z1 = c*h4;			t3 = *zki;	pzki1++;
	    z2 = s*h4;			zki++;		pzki++;
	    z3 = s*t4;	s1 += s3;	h4 = *(zki1-1);	*pzki1 = q1;	
	    z4 = c*t4;	s4 -= s2;			*pzki = q4;  	
	    } while(k > 4);
	    // } while(k > 0);
	//
	if(k > 0)
	    {
	    r1 = c*h1;			t4 = *zki;	k -= 4;
	    r2 = s*h1;			zki++;       	
	    r3 = s*t1;	z1 += z3;		        *(pzki1+1) = s1;	
	    r4 = c*t1;	z4 -= z2;	zki1++;      	*(pzki+1) = s4;     
	    q1 = c*h2;			      
	    q2 = s*h2;			zki++;       	
	    q3 = s*t2;	r1 += r3;			*(pzki1+2) = z1;	
	    q4 = c*t2;	r4 -= r2;	zki1++;      	*(pzki+2) = z4;  	
	    s1 = c*h3;					pzki1 = zki1 - 6;
	    s2 = s*h3;			zki++; 		pzki = zki - 6;
	    s3 = s*t3;	q1 += q3;			*pzki1 = r1;	
	    s4 = c*t3;	q4 -= q2;	zki1+= 2;	*pzki = r4;
	    z1 = c*h4;					pzki1++;
	    z2 = s*h4;			zki++;		pzki++;
	    z3 = s*t4;	s1 += s3;			*pzki1 = q1;	
	    z4 = c*t4;	s4 -= s2;			*pzki = q4;  	
	    }
	//

							pzki1++;     
							pzki++;      
			z1 += z3;			*pzki1 = s1;	
			z4 -= z2;			*pzki = s4;     
							pzki1++;     
							pzki++;      
							*pzki1 = z1;	
							*pzki = z4;  	
	}
    zki1 = pzki1+1;
    zki = pzki+1;
    for(k = 0; k < nr; k++)
	{
	h1 = *zki1;	t1 = *zki;
	r1 = c*h1;	pzki1 = zki1;
	r2 = s*h1;	pzki = zki;
	r3 = s*t1;	zki++; 
	r4 = c*t1;	zki1++;
	r1 += r3;
	r4 -= r2;
	*pzki1 = r1;
	*pzki = r4;
	}
    }

static double extvar = 0.0;

double dotprod(double *a,double *b, int l)
    {
    register double s0 = 0, s1 = 0, s2 = 0, s3 = 0;	// 4
    register double t0 = 0, t1 = 0, t2 = 0, t3 = 0;	// 4
    register double a0,b0,a1,b1,a2,b2,a3,b3;
    register double prea,preb;
    int i,ll = l&3;
    for(i = l-ll; i < l; i++)
	s3 += a[i] * b[i];
    l -= ll;
    int lend = l - 8;
    if(lend < 0) lend = 0;
    for(i = 0; i < lend; i += 4)
	{
	s0 += t0;	a0 = a[0];	b0 = b[0];
	s1 += t1;	a1 = a[1];	b1 = b[1];
	s2 += t2;	a2 = a[2];	b2 = b[2];
	s3 += t3;	a3 = a[3];	b3 = b[3];
	t0 = a0 * b0;	a += 4;		b += 4;
	t1 = a1 * b1;	prea = a[4];	preb = b[4];
	t2 = a2 * b2;
	t3 = a3 * b3;
	}
    for(i = lend; i < l; i += 4)
	{
	s0 += t0;	a0 = a[0];	b0 = b[0];
	s1 += t1;	a1 = a[1];	b1 = b[1];
	s2 += t2;	a2 = a[2];	b2 = b[2];
	s3 += t3;	a3 = a[3];	b3 = b[3];
	t0 = a0 * b0;	a += 4;		b += 4;
	t1 = a1 * b1;
	t2 = a2 * b2;
	t3 = a3 * b3;
	}
    s0 += t0;
    s1 += t1;
    s2 += t2;
    s3 += t3;
    if(l == -1)
	extvar = prea + preb;
    return s0+s1+s2+s3;
    }


double dotprod2(double *a,int inca, double *b, int incb, int l)
    {
    register double s0 = 0, s1 = 0, s2 = 0, s3 = 0;	// 4
    register double t0 = 0, t1 = 0, t2 = 0, t3 = 0;	// 4
    register double a0,b0,a1,b1,a2,b2,a3,b3;
    register double *pa0,*pa1,*pa2,*pa3, *pb0,*pb1,*pb2,*pb3;
    int i,ll = l&3;
    int lmain = l - ll;
    if(lmain < 0) lmain = 0;
    pa0 = a;
    pa1 = pa0 + inca;
    pa2 = pa1 + inca;
    pa3 = pa2 + inca;
    pb0 = b;
    pb1 = pb0 + incb;
    pb2 = pb1 + incb;
    pb3 = pb2 + incb;
    int inca4 = inca << 2;
    int incb4 = incb << 2;
    i = lmain;
    while(i > 0)
	{
	s0 += t0;	a0 = *pa0;	b0 = *pb0;
	s1 += t1;	a1 = *pa1;	b1 = *pb1;
	i -= 4;
	s2 += t2;	a2 = *pa2;	b2 = *pb2;
	s3 += t3;	a3 = *pa3;	b3 = *pb3;
	t0 = a0 * b0;	pa0 += inca4;	pb0 += incb4;
	t1 = a1 * b1;	pa1 += inca4;	pb1 += incb4;
	t2 = a2 * b2;	pa2 += inca4;	pb2 += incb4;
	t3 = a3 * b3;	pa3 += inca4;	pb3 += incb4;
	} 
    s0 += t0;
    s1 += t1;
    s2 += t2;
    s3 += t3;
    for(i = lmain; i < l; i++)
	{
	a0 = *pa0;	b0 = *pb0;
	pa0 += inca;	pb0 += incb;
	s0 += a0 * b0;	
	}
    return s0+s1+s2+s3;
    }

// A = U d V; Columns of U are orthog; Rows of V are; d is diagonal
// Dimensions are  (m x n) = (m x k) (k x k) (k x n) where k = min(m,n)
//
/*
void SVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V)
    {
    Matrix rho;
    if(A.Nrows() <= A.Ncols())
	rho = A * A.t();
    else
	rho = A.t() * A;
    Vector evals;
    rho *= -1.0;
    if(A.Nrows() <= A.Ncols())
	{
	EigenValues(rho,evals,U);
	V = U.t() * A;
	d.ReDimension(A.Nrows());
	for(int i = 1; i <= d.Length(); i++)
	    {
	    d(i) = Norm(V.Row(i));
	    if(d(i) != 0.0) V.Row(i) *= 1.0 / d(i);
	    }
	}
    else
	{
	Matrix Vt;
	EigenValues(rho,evals,Vt);
	V = Vt.t();
	U = A * Vt;
	d.ReDimension(A.Ncols());
	for(int i = 1; i <= d.Length(); i++)
	    {
	    d(i) = Norm(U.Column(i));
	    if(d(i) != 0.0) U.Column(i) *= 1.0 / d(i);
	    }
	}
    }
*/

// A = U d V; Columns of U are orthog; Rows of V are; d is diagonal
// Dimensions are  (m x n) = (m x k) (k x k) (k x n) where k = min(m,n)

void dosvd(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V);

/*
void SVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V)
    {
    if(A.Nrows() < A.Ncols())
	{
	Matrix AA = A.t();
	Matrix Ut, Vt;
	SVD(AA,Ut,d,Vt);
	U = Vt.t();
	V = Ut.t();
	return;
	}
    U.ReDimension(A.Nrows(),A.Ncols());
    d.ReDimension(A.Ncols());
    V.ReDimension(A.Ncols(),A.Ncols());
    Matrix UU(U),VV;
    UU.ReDimension(A.Nrows(),A.Ncols());
    VV.ReDimension(A.Ncols(),A.Ncols());
    dosvd(A,UU,d,VV);
    Vector dd(-d);
    IntArray1 ind(1,d.Length());
    int i;
    for(i = 1; i <= d.Length(); i++)
	ind[i] = i;
    Sort(dd,ind);
    d = -dd;
    for(i = 1; i <= UU.Ncols(); i++)
	U.Column(i) = UU.Column(ind(i));
    for(i = 1; i <= VV.Nrows(); i++)
	V.Row(i) = VV.Row(ind(i));
    }
*/

int svd(double *a, int *m, int *n, int *mp, int *np, 
	double *w, double *v);

void dosvd(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V)
    {
    Matrix AA = A.t();
    int m = A.Nrows();
    int n = A.Ncols();
    Real maxel = 0.0;
    for(int i = 1; i <= m; i++)
	for(int j = 1; j <= n; j++)
	    if(fabs(AA(j,i)) > maxel)
		maxel = fabs(AA(j,i));

    if(maxel != 0.0) AA *= 1.0 / maxel;
    svd(AA.Store(),&m,&n,&m,&n,d.Store(),V.Store());
    if(maxel != 0.0) d *= maxel;
    U = AA.t();
    }

/* svd.f -- translated by f2c (version 19951025).
*/

// #include "f2c.h"

/* Table of constant values */

static double one = 1.;

inline double d_sign(double *x, double *y)
    {
    return (*y >= 0 ? fabs(*x) : -fabs(*x));
    }

double dpythag_(double *a, double *b);

int svd(double *a, int *m, int *n, int *mp, int *np, 
	double *w, double *v)
    {
    int a_dim1, a_offset, v_dim1, v_offset, i__1, i__2, i__3;
    double d__1, d__2, d__3, d__4;

    static double c__, f, g, h__;
    static int i__, j, k, l;
    static double s, scale, x, y, z__, anorm;
    static int jj, nm;
    static double rv1[2000];
    static int its;

    v_dim1 = *np;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    --w;
    a_dim1 = *mp;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    g = 0.;
    scale = 0.;
    anorm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) 
	{
	l = i__ + 1;
	rv1[i__ - 1] = scale * g;
	g = 0.;
	s = 0.;
	scale = 0.;
	if (i__ <= *m) 
	    {
	    i__2 = *m;
	    for (k = i__; k <= i__2; ++k) 
		scale += (d__1 = a[k + i__ * a_dim1], fabs(d__1));
	    if (scale != 0.) 
		{
		i__2 = *m;
		for (k = i__; k <= i__2; ++k) 
		    {
		    a[k + i__ * a_dim1] /= scale;
		    s += a[k + i__ * a_dim1] * a[k + i__ * a_dim1];
		    }
		f = a[i__ + i__ * a_dim1];
		d__1 = sqrt(s);
		g = -d_sign(&d__1, &f);
		h__ = f * g - s;
		a[i__ + i__ * a_dim1] = f - g;
		i__2 = *n;
		for (j = l; j <= i__2; ++j) 
		    {
		    s = 0.;
		    i__3 = *m;
		    for (k = i__; k <= i__3; ++k) 
			s += a[k + i__ * a_dim1] * a[k + j * a_dim1];
		    f = s / h__;
		    i__3 = *m;
		    for (k = i__; k <= i__3; ++k) 
			a[k + j * a_dim1] += f * a[k + i__ * a_dim1];
		    }
		i__2 = *m;
		for (k = i__; k <= i__2; ++k) 
		    a[k + i__ * a_dim1] = scale * a[k + i__ * a_dim1];
		}
	    }
	w[i__] = scale * g;
	g = 0.;
	s = 0.;
	scale = 0.;
	if (i__ <= *m && i__ != *n) 
	    {
	    i__2 = *n;
	    for (k = l; k <= i__2; ++k) 
		scale += (d__1 = a[i__ + k * a_dim1], fabs(d__1));
	    if (scale != 0.) 
		{
		i__2 = *n;
		for (k = l; k <= i__2; ++k) 
		    {
		    a[i__ + k * a_dim1] /= scale;
		    s += a[i__ + k * a_dim1] * a[i__ + k * a_dim1];
		    }
		f = a[i__ + l * a_dim1];
		d__1 = sqrt(s);
		g = -d_sign(&d__1, &f);
		h__ = f * g - s;
		a[i__ + l * a_dim1] = f - g;
		i__2 = *n;
		for (k = l; k <= i__2; ++k) 
		    rv1[k - 1] = a[i__ + k * a_dim1] / h__;
		i__2 = *m;
		for (j = l; j <= i__2; ++j) 
		    {
		    s = 0.;
		    i__3 = *n;
		    for (k = l; k <= i__3; ++k) 
			s += a[j + k * a_dim1] * a[i__ + k * a_dim1];
		    i__3 = *n;
		    for (k = l; k <= i__3; ++k) 
			a[j + k * a_dim1] += s * rv1[k - 1];
		    }
		i__2 = *n;
		for (k = l; k <= i__2; ++k) 
		    a[i__ + k * a_dim1] = scale * a[i__ + k * a_dim1];
		}
	    }
	// Computing MAX
	d__3 = anorm, d__4 = (d__1 = w[i__], fabs(d__1)) + (d__2 = rv1[i__ - 1]
		, fabs(d__2));
	anorm = max(d__3,d__4);
	}
    for (i__ = *n; i__ >= 1; --i__) 
	{
	if (i__ < *n) 
	    {
	    if (g != 0.) 
		{
		i__1 = *n;
		for (j = l; j <= i__1; ++j) 
		    v[j + i__ * v_dim1] = a[i__ + j * a_dim1] / a[i__ + l * 
				a_dim1] / g;
		i__1 = *n;
		for (j = l; j <= i__1; ++j) 
		    {
		    s = 0.;
		    i__2 = *n;
		    for (k = l; k <= i__2; ++k) 
			s += a[i__ + k * a_dim1] * v[k + j * v_dim1];
		    i__2 = *n;
		    for (k = l; k <= i__2; ++k) 
			v[k + j * v_dim1] += s * v[k + i__ * v_dim1];
		    }
		}
	    i__1 = *n;
	    for (j = l; j <= i__1; ++j) 
		{
		v[i__ + j * v_dim1] = 0.;
		v[j + i__ * v_dim1] = 0.;
		}
	    }
	v[i__ + i__ * v_dim1] = 1.;
	g = rv1[i__ - 1];
	l = i__;
	}
    for (i__ = min(*m,*n); i__ >= 1; --i__) 
	{
	l = i__ + 1;
	g = w[i__];
	i__1 = *n;
	for (j = l; j <= i__1; ++j) 
	    a[i__ + j * a_dim1] = 0.;
	if (g != 0.) 
	    {
	    g = 1. / g;
	    i__1 = *n;
	    for (j = l; j <= i__1; ++j) 
		{
		s = 0.;
		i__2 = *m;
		for (k = l; k <= i__2; ++k) 
		    s += a[k + i__ * a_dim1] * a[k + j * a_dim1];
		f = s / a[i__ + i__ * a_dim1] * g;
		i__2 = *m;
		for (k = i__; k <= i__2; ++k) 
		    a[k + j * a_dim1] += f * a[k + i__ * a_dim1];
		}
	    i__1 = *m;
	    for (j = i__; j <= i__1; ++j) 
		a[j + i__ * a_dim1] *= g;
	    } 
	else 
	    {
	    i__1 = *m;
	    for (j = i__; j <= i__1; ++j) 
		a[j + i__ * a_dim1] = 0.;
	    }
	a[i__ + i__ * a_dim1] += 1.;
	}
    for (k = *n; k >= 1; --k) 
	{
	for (its = 1; its <= 30; ++its) 
	    {
	    for (l = k; l >= 1; --l) 
		{
		nm = l - 1;
		if ((d__1 = rv1[l - 1], fabs(d__1)) + anorm == anorm) 
		    goto L2;
		if ((d__1 = w[nm], fabs(d__1)) + anorm == anorm) 
		    goto L1;
		}
L1:
	    c__ = 0.;
	    s = 1.;
	    i__1 = k;
	    for (i__ = l; i__ <= i__1; ++i__) 
		{
		f = s * rv1[i__ - 1];
		rv1[i__ - 1] = c__ * rv1[i__ - 1];
		if (fabs(f) + anorm == anorm) 
		    {
		    goto L2;
		    }
		g = w[i__];
		h__ = dpythag_(&f, &g);
		w[i__] = h__;
		h__ = 1. / h__;
		c__ = g * h__;
		s = -(f * h__);
		i__2 = *m;
		for (j = 1; j <= i__2; ++j) 
		    {
		    y = a[j + nm * a_dim1];
		    z__ = a[j + i__ * a_dim1];
		    a[j + nm * a_dim1] = y * c__ + z__ * s;
		    a[j + i__ * a_dim1] = -(y * s) + z__ * c__;
		    }
		}
L2:
	    z__ = w[k];
	    if (l == k) 
		{
		if (z__ < 0.) 
		    {
		    w[k] = -z__;
		    i__1 = *n;
		    for (j = 1; j <= i__1; ++j) 
			v[j + k * v_dim1] = -v[j + k * v_dim1];
		    }
		goto L3;
		}
	    if (its >= 30) 
		{
		cout << "No convergence in svd\n";
		}
	    x = w[l];
	    nm = k - 1;
	    y = w[nm];
	    g = rv1[nm - 1];
	    h__ = rv1[k - 1];
	    f = ((y - z__) * (y + z__) + (g - h__) * (g + h__)) / (h__ * 2. * 
		    y);
	    g = dpythag_(&f, &one);
	    f = ((x - z__) * (x + z__) + h__ * (y / (f + d_sign(&g, &f)) - 
		    h__)) / x;
	    c__ = 1.;
	    s = 1.;
	    i__1 = nm;
	    for (j = l; j <= i__1; ++j) 
		{
		i__ = j + 1;
		g = rv1[i__ - 1];
		y = w[i__];
		h__ = s * g;
		g = c__ * g;
		z__ = dpythag_(&f, &h__);
		rv1[j - 1] = z__;
		c__ = f / z__;
		s = h__ / z__;
		f = x * c__ + g * s;
		g = -(x * s) + g * c__;
		h__ = y * s;
		y *= c__;
		i__2 = *n;
		for (jj = 1; jj <= i__2; ++jj) 
		    {
		    x = v[jj + j * v_dim1];
		    z__ = v[jj + i__ * v_dim1];
		    v[jj + j * v_dim1] = x * c__ + z__ * s;
		    v[jj + i__ * v_dim1] = -(x * s) + z__ * c__;
		    }
		z__ = dpythag_(&f, &h__);
		w[j] = z__;
		if (z__ != 0.) 
		    {
		    z__ = 1. / z__;
		    c__ = f * z__;
		    s = h__ * z__;
		    }
		f = c__ * g + s * y;
		x = -(s * g) + c__ * y;
		i__2 = *m;
		for (jj = 1; jj <= i__2; ++jj) 
		    {
		    y = a[jj + j * a_dim1];
		    z__ = a[jj + i__ * a_dim1];
		    a[jj + j * a_dim1] = y * c__ + z__ * s;
		    a[jj + i__ * a_dim1] = -(y * s) + z__ * c__;
		    }
		}
	    rv1[l - 1] = 0.;
	    rv1[k - 1] = f;
	    w[k] = x;
	    }
L3:
	;
	}
    return 0;
    }

double dpythag_(double *a, double *b)
    {
    double ret_val, d__1;
    static double absa, absb;

    absa = fabs(*a);
    absb = fabs(*b);
    if (absa > absb) 

	{
	d__1 = absb / absa;
	ret_val = absa * sqrt(d__1 * d__1 + 1.);
	} 
    else 
	{
	if (absb == 0.) 
	    ret_val = 0.;
	else 
	    {
	    d__1 = absa / absb;
	    ret_val = absb * sqrt(d__1 * d__1 + 1.);
	    }
	}
    return ret_val;
    } 


typedef long int lint;


void getrowbasis(Matrix& B, Real cutoff,bool nopivot = false)	// cutoff is in norm^2 of a row, 10^-20 is OK
    {
    int i, m = B.Nrows();
    Vector v(B.Ncols()), orignorm(m);
    for(i = 1; i <= m; i++)
	orignorm(i) = B.Row(i) * B.Row(i);
    for(i = 1; i <= m; i++)
	{
	int k = 1;
	Real maxnor = 0.0,nor, sw;
	if(nopivot)
	    k = i, maxnor = B.Row(i) * B.Row(i);
	else
	    for(int j = i; j <= m; j++)
		if((nor = B.Row(j)*B.Row(j)) > maxnor)
		    k = j, maxnor = nor;
	if(maxnor <= cutoff) break;
	if(maxnor < orignorm(k)*1.0e-4)	// probable loss of orthogonality
	    {
	    for(int j = 1; j < i; j++)
		B.Row(k) += (-1.0 * (B.Row(k) * B.Row(j))) * B.Row(j);
	    maxnor = B.Row(k) * B.Row(k);
	    }
	B.Row(k) *= 1.0 / sqrt(maxnor);
	v = B.Row(k);
	if(k != i) 
	    {
	    B.Row(k) = B.Row(i); B.Row(i) = v; sw = orignorm(i); 
	    orignorm(i) = orignorm(k), orignorm(k)=sw;
	    }
	for(int j = i+1; j <= m; j++)
	    B.Row(j) += ((-1.0) * (v * B.Row(j))) * v;
	}
    if(--i == 0) i = 1, B(1,1) = 1.0;
    B.ReduceDimension(i,B.Ncols());
    }

void startSVD(const Matrix& A,Matrix& R, Matrix &U,Real cutsq = 1.0e-20,bool nopivot = false)
    {
    U = A;
    getrowbasis(U,cutsq);
    R = A * U.t();
    //cout << "R is\n" << R;
    }

void onepassSVD(const Matrix& A,Matrix& U, Matrix& D, Matrix& V1, Matrix &V2,Real cutsq = 1.0e-20)
    {
    Matrix evecs,r1,r2;
    Vector evals;
    startSVD(A,r1,V2,cutsq);
    Matrix rsq = (-1.0) * r1 * r1.t();
    EigenValues(rsq,evals,evecs);
    U = evecs;
    r2 = evecs.t() * r1;
    startSVD(r2,D,V1);
    }

/*
void thinsvd(const Matrix& A,Matrix& U, Vector& d,  Matrix &V,Real cutsq)
    {
    int m = A.Nrows(), n = A.Ncols();
    //cout << "starting thinsvd, m, n are " << m SP n << endl;
    if(m > n)
	{
	Matrix At = A.t(), Ut, Vt;
	thinsvd(At,Ut,d,Vt,cutsq);
	U = Vt.t();
	V = Ut.t();
	return;
	}
    Matrix VV, r1;
    startSVD(A,r1,VV,cutsq,false);
    int k = r1.Ncols();
    if(k == 1)
	{
	d.ReDimension(1);
	d(1) = Norm(r1.TreatAsVector());
	U = r1;
	if(d(1) == 0.0)
	    U = 0.0, U(1,1) = 1.0;
	else
	    U *= 1.0 / d(1);
	V = VV;
	return;
	}
    if(k < m)
	{
	//cout << "k, m are " << k SP m << ", calling recursively" << endl;
	Matrix r1t = r1.t(), Ut, Vt;
	thinsvd(r1t,Ut,d,Vt,-1.0);
	U = Vt.t();
	V = Ut.t() * VV;
	return;
	}
    Matrix vv;
    newSVD(r1,U,d,vv);
    V = vv * VV;
    }
    */

void checkSVD(const Matrix& A,const Matrix &U, const Vector &d, const Matrix& V)
    {
    int m = d.Length();
    Matrix D(m,m); D = 0.0; D.Diagonal() = d;
    Matrix AA = U * D * V;
    Matrix checkit = AA - A;
    Real err = Norm(checkit.TreatAsVector());
    err *= 1.0 / Norm(A.TreatAsVector());
    Matrix chu = U.t() * U;
    chu += -1.0;
    Real erru = Norm(chu.TreatAsVector());
    erru *= 1.0 / (chu.Nrows() * chu.Ncols());
    Matrix chv = V * V.t();
    chv += -1.0;
    Real errv = Norm(chv.TreatAsVector());
    errv *= 1.0 / (chv.Nrows() * chv.Ncols());
    cout << "ave errors in SVD * 10^14 is " << err * 1.0e14 SP erru * 1.0e14 SP errv * 1.0e14 << endl;
    }
Real check_complex_SVD(const Matrix& Mre, const Matrix& Mim, const Matrix& Ure, 
	        const Matrix& Uim, const Vector& d, const Matrix& Vre, const Matrix& Vim)
    {
    Matrix ddi(d.Length(),d.Length());
    ddi = 0.0;
    ddi.Diagonal() = d;
    Matrix resre = Ure * ddi * Vre - Uim * ddi * Vim;
    Matrix resim = Ure * ddi * Vim + Uim * ddi * Vre;
    return Norm(Matrix(Mre-resre).TreatAsVector()) + Norm(Matrix(Mim-resim).TreatAsVector());
    }

void newcomplexSVD(const Matrix& Are, const Matrix& Aim, Matrix& Ure, Matrix& Uim, 
	Vector& d, Matrix& Vre, Matrix& Vim)
    {
    LAPACK_INT m = Are.Nrows(), n = Are.Ncols(); 
    if(m < n)
	{
	Matrix mret = Are.t(), mimt = -Aim.t(),UUre,UUim,VVre,VVim;
	newcomplexSVD (mret,mimt, UUre, UUim, d, VVre, VVim);
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
	    AA(i,2*j-1) = Are(j,i), AA(i,2*j) = Aim(j,i);

    Matrix UU(n,2*m), VV(n,2*n);
    Ure.ReDimension(m,n);
    Uim.ReDimension(m,n);
    Vre.ReDimension(n,n);
    Vim.ReDimension(n,n);
    d.ReDimension(n);
    LAPACK_INT lda = m, ldu = m, ldv = n, info;
    LAPACK_INT lwork = -1;
    Vector work(2*n*n+2*n+m+100), iwork(8*n), rwork(10*n*n+10 *n+1000);
    lwork = n*n+n+m+50;

    /*
    zgesdd_(&jobz,&m,&n,(LAPACK_COMPLEX*)AA.Store(), 
	    &lda,d.Store(), (LAPACK_COMPLEX*)UU.Store(), 
	    &ldu, (LAPACK_COMPLEX*)VV.Store(), &ldv, 
	    (LAPACK_COMPLEX*)work.Store(), &lwork, (LAPACK_INT *)iwork.Store(), &info);
     JOBZ, [M], [N], A, [LDA], S, U, [LDU], VT, [LDVT], 
		  *       [WORK], [LWORK], [RWORK], [IWORK], [INFO])
    int zgesdd_(char*, LAPACK_INT*, LAPACK_INT*, LAPACK_COMPLEX*,
	    LAPACK_INT*, LAPACK_REAL*, LAPACK_COMPLEX*,
	    LAPACK_INT*, LAPACK_COMPLEX*, LAPACK_INT*,
	    LAPACK_COMPLEX*, LAPACK_INT*, LAPACK_REAL*,
	    LAPACK_INT*, LAPACK_INT*)' 
    lwork = (int)work(1);
    //cout << "optimal size is " << work(1) << endl;
    //cout << "m, n, lwork are " << m SP n SP lwork << endl;
    work.ReDimension(lwork*2);
    */

    zgesdd_wrapper(&jobz,&m,&n,(LAPACK_COMPLEX*)AA.Store(), &lda,d.Store(), (LAPACK_COMPLEX*)UU.Store(), &ldu,
	    (LAPACK_COMPLEX*)VV.Store(), &ldv, (LAPACK_COMPLEX*)work.Store(), &lwork,
	    rwork.Store(), (LAPACK_INT *)iwork.Store(), &info);
    bool do_over = false;
    if(info != 0) do_over = true;
    //dgesdd_(&jobz,&m,&n,AA.Store(), &lda,d.Store(), UU.Store(), &ldu,
//	    VV.Store(), &ldv, work.Store(), &lwork,
//	    (LAPACK_INT *)iwork.Store(), &info);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
	    Ure(j,i) = UU(i,2*j-1), Uim(j,i) = UU(i,2*j);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= n; j++)
	    Vre(j,i) = VV(i,2*j-1), Vim(j,i) = VV(i,2*j);
    if(!do_over)
	{
	Real err = check_complex_SVD( Are,  Aim,  Ure, Uim,  d,  Vre,  Vim);
	if(err > (Norm(Are.TreatAsVector()) + Norm(Aim.TreatAsVector())) * 1.0e-8)
	    {
	    do_over = true;
	    cout << "doing newcomplexSVD over, zgesdd err is " << err << endl;
	    }
	}
    if(!do_over) return;
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
	    AA(i,2*j-1) = Are(j,i), AA(i,2*j) = Aim(j,i);
    Matrix r(AA);
    r.Randomize();
    AA += (Norm(AA.TreatAsVector()) * 1.0e-10) * r; 
    zgesvd_(&jobz,&jobz,&m,&n,(LAPACK_COMPLEX*)AA.Store(), &lda,d.Store(), (LAPACK_COMPLEX*)UU.Store(), &ldu,
	    (LAPACK_COMPLEX*)VV.Store(), &ldv, (LAPACK_COMPLEX*)work.Store(), &lwork,
	    rwork.Store(), &info);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
	    Ure(j,i) = UU(i,2*j-1), Uim(j,i) = UU(i,2*j);
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= n; j++)
	    Vre(j,i) = VV(i,2*j-1), Vim(j,i) = VV(i,2*j);
    do_over = false;
    if(info != 0) do_over = true;
    if(!do_over)
	{
	Real err = check_complex_SVD( Are,  Aim,  Ure, Uim,  d,  Vre,  Vim);
	if(err > (Norm(Are.TreatAsVector()) + Norm(Aim.TreatAsVector())) * 1.0e-5)
	    {
	    cout << "messed up newcomplexSVD, zgesvd err is " << err << endl;
	    do_over = true;
	    }
	}
    if(!do_over) return;
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= m; j++)
	    AA(i,2*j-1) = Are(j,i), AA(i,2*j) = Aim(j,i);
    Matrix re(Are), im(Aim), Mre(Are), Mim(Aim);
    re.Randomize();
    im.Randomize();
    Mre += (Norm(Are.TreatAsVector()) * 1.0e-11) * re; 
    Mim += (Norm(Aim.TreatAsVector()) * 1.0e-11) * im; 
    void SVDcomplex(const Matrix& Mre, const Matrix& Mim, Matrix& Ure, 
	        Matrix& Uim, Vector& d, Matrix& Vre, Matrix& Vim);
    cout << "Norm of Mre is " << Norm(Mre.TreatAsVector()) << endl;
    cout << "Norm of Mim is " << Norm(Mim.TreatAsVector()) << endl;
    SVDcomplex(Mre,Mim,Ure,Uim,d,Vre,Vim);

    if(1)
	{
	Real err = check_complex_SVD( Are,  Aim,  Ure, Uim,  d,  Vre,  Vim);
	if(err > (Norm(Are.TreatAsVector()) + Norm(Aim.TreatAsVector())) * 1.0e-5)
	    {
	    cout << "messed up newcomplexSVD SVDcomplex, err is " << err << endl;
	    cout << "d is " << d << endl;
	    cout << "Norm of Are is " << Norm(Are.TreatAsVector()) << endl;
	    cout << "Norm of Aim is " << Norm(Aim.TreatAsVector()) << endl;
	    cout << "Norm of Ure is " << Norm(Ure.TreatAsVector()) << endl;
	    cout << "Norm of Uim is " << Norm(Uim.TreatAsVector()) << endl;
	    cout << "Norm of Vre is " << Norm(Vre.TreatAsVector()) << endl;
	    cout << "Norm of Vim is " << Norm(Vim.TreatAsVector()) << endl;
	    std::ofstream outfile;
	    outfile.open("bad_matrix");
	    int nrows = Are.Nrows(), ncols = Are.Ncols();
	    outfile.write((char*)&nrows,sizeof(nrows)); 
	    outfile.write((char*)&ncols,sizeof(ncols)); 
	    outfile.write((char*)Are.Store(),sizeof(Real)*nrows*ncols); 
	    outfile.write((char*)Aim.Store(),sizeof(Real)*nrows*ncols); 
	    nrows = Ure.Nrows(), ncols = Ure.Ncols();
	    outfile.write((char*)&nrows,sizeof(nrows)); 
	    outfile.write((char*)&ncols,sizeof(ncols)); 
	    outfile.write((char*)Ure.Store(),sizeof(Real)*nrows*ncols); 
	    outfile.write((char*)Uim.Store(),sizeof(Real)*nrows*ncols); 
	    nrows = Vre.Nrows(), ncols = Vre.Ncols();
	    outfile.write((char*)&nrows,sizeof(nrows)); 
	    outfile.write((char*)&ncols,sizeof(ncols)); 
	    outfile.write((char*)Vre.Store(),sizeof(Real)*nrows*ncols); 
	    outfile.write((char*)Vim.Store(),sizeof(Real)*nrows*ncols); 
	    nrows = d.Length();
	    outfile.write((char*)&nrows,sizeof(nrows)); 
	    outfile.write((char*)d.Store(),sizeof(Real)*nrows); 
	    outfile.close();
	    error("bad newcomplexSVD");
	    }
	}
    }


inline double abssq(Complex a)
    {
    return real(conj(a) * a);
    }

void CSVD(ComplexMatrix a,   ComplexMatrix& u, Vector& s, ComplexMatrix& v)
/* Singular Value Decomposition, a = u * s * Conj(Tran(v)), a is destroyed by CSVD
   the diagonal matrix s is output as a vector, m must be >= n
   if smaller, a should be filled with zero rows
   this code is adapted from Collected Algorithms from ACM, Algorithm 358
   The transformation Conj(Tran(u)) is applied to the p vectors given in columns
   n, n+1, ..., n+p-1 of matrix a
   See: http://www.scs.fsu.edu/~burkardt/f77_src/toms358/toms358.f
   and: http://www.scs.fsu.edu/~burkardt/f77_src/toms358/toms358_prb.f
*/
    {
    int m = a.nrow, n = a.ncol, p = 0, nu = n, nv = n;
    cout << "in CSVD, m,n are " << m SP n << endl;
    s.ReDimension(n);
    u = ComplexMatrix(m,n);
    v = ComplexMatrix(n,n);
    int maxA = max(m,n) + 100;
    Vector b(maxA), c(maxA), t(maxA);
    double cs, eps, eta, f, g, h;
    int i, j, k, k1, L, L1, nM1, np;
    Complex q, r;
    double sn;
    double tol, w, x, y, z;
    eta = 2.8E-16;			/* eta = the relative machine precision */
    tol = 4.0E-293; 		/* tol = the smallest normalized positive number, divided by eta */
    /* eta = 2^-52 * 1.26 fudge*/
    /* tol = 2^-1023 / eta */
    np = n + p;
    nM1 = n - 1;
    L = 0;
    /*
       HOUSEHOLDER REDUCTION
     */
    c.el(0) = 0.0;
    k = 0;
    while (1)
	{
	cout << "place 1, k = " << k << endl;
	k1 = k + 1;
	/*
	   ELIMINATION OF a.el(i,k), i = k, ..., m-1
	 */
	z = 0.0;
	for (i = k; i < m; i++)
	    z += abssq(a.el(i,k));
	b.el(k) = 0.0;
	if (z > tol)
	    {
	    z = sqrt(z);
	    b.el(k) = z;
	    w = abs(a.el(k,k));
	    q = 1.0;
	    if (w != 0.0) q = a.el(k,k) / w;
	    a.el(k,k) = q * (z + w);
	    if (k != np - 1)
		{
		for (j = k1; j < np; j++)
		    {
		    q = 0.0;
		    for (i = k; i < m; i++)
			q += conj(a.el(i,k)) * a.el(i,j);
		    q /= z * (z + w);
		    for (i = k; i < m; i++)
			a.el(i,j) -= q * a.el(i,k);
		    }
		}
	    /*
	       PHASE TRANSFORMATION
	     */
	    q = -conj(a.el(k,k)) / abs(a.el(k,k));
	    for (j = k1; j < np; j++)
		a.el(k,j) *= q;
	    }
	/*
	   ELIMINATION OF a.el(k,j), j = k+2, ..., n-1
	 */
	cout << "place 2, k, nM1 = " << k SP nM1 << endl;
	if (k == nM1) break;
	z = 0.0;
	for (j = k1; j < n; j++)
	    z += abssq(a.el(k,j));
	c.el(k1) = 0.0;
	if (z > tol)
	    {
	    z = sqrt(z);
	    c.el(k1) = z;
	    w = abs(a.el(k,k1));
	    q = 1.0;
	    if (w != 0.0) q = a.el(k,k1) / w;
	    a.el(k,k1) = q * (z + w);
	    for (i = k1; i < m; i++)
		{
		q = 0.0;
		for (j = k1; j < n; j++)
		    q += conj(a.el(k,j)) * a.el(i,j);
		q /= z * (z + w);
		for (j = k1; j < n; j++) 
		    a.el(i,j) -= q * a.el(k,j);
		}
	    /*
	       PHASE TRANSFORMATION
	     */
	    q = -conj(a.el(k,k1)) / abs(a.el(k,k1));
	    for (i = k1; i < m; i++)
		a.el(i,k1) *= q;
	    }
	k = k1;
	}
    /*
       TOLERANCE FOR NEGLIGIBLE ELEMENTS
     */
    eps = 0.0;
    cout << "place 3 " <<endl;
    for (k = 0; k < n; k++)
	{
	s.el(k) = b.el(k);
	t.el(k) = c.el(k);
	if (s.el(k) + t.el(k) > eps)
	    eps = s.el(k) + t.el(k);
	}
    eps *= eta;
    /*
       INITIALIZATION OF u AND v
     */
    for (j = 0; j < nu; j++)
	{
	for (i = 0; i < m; i++)
	    u.el(i,j) = 0.0;
	u.el(j,j) = 1.0;
	}
    for (j = 0; j < nv; j++)
	{
	for (i = 0; i < n; i++)
	    v.el(i,j) = 0.0;
	v.el(j,j) = 1.0;
	}
    /*
       QR DIAGONALIZATION
     */
    for (k = nM1; k >= 0; k--)
	{
	cout << "place 4 , k = " << k <<endl;
	/*
	   TEST FOR SPLIT
	 */
	while (1)
	    {
	    for (L = k; L >= 0; L--)
		{
		if (abs(t.el(L)) <= eps) goto Test;
		if (abs(s.el(L - 1)) <= eps) break;
		}
	    /*
	       CANCELLATION OF E(L)
	     */
	    cs = 0.0;
	    sn = 1.0;
	    L1 = L - 1;
	    for (i = L; i <= k; i++)
		{
		f = sn * t.el(i);
		t.el(i) *= cs;
		if (abs(f) <= eps) goto Test;
		h = s.el(i);
		w = sqrt(f * f + h * h);
		s.el(i) = w;
		cs = h / w;
		sn = -f / w;
		for (j = 0; j < n; j++)
		    {
		    x = real(u.el(j,L1));
		    y = real(u.el(j,i));
		    u.el(j,L1) = x * cs + y * sn;
		    u.el(j,i) = y * cs - x * sn;
		    }
		}
	    /*
	       TEST FOR CONVERGENCE
	     */
Test:	    w = s.el(k);
	cout << "place 5 , k, L = " << k SP L <<endl;
	    if (L == k) break;
	/*
	   ORIGIN SHIFT
	 */
	    x = s.el(L);
	    y = s.el(k - 1);
	    g = t.el(k - 1);
	    h = t.el(k);
	    f = ((y - w) * (y + w) + (g - h) * (g + h)) / (2.0 * h * y);
	    g = sqrt(f * f + 1.0);
	    if (f < 0.0) g = -g;
	    f = ((x - w) * (x + w) + (y / (f + g) - h) * h) / x;
	    /*
	       QR STEP
	     */
	    cs = 1.0;
	    sn = 1.0;
	    L1 = L + 1;
	    for (i = L1; i <= k; i++)
		{
		g = t.el(i);
		y = s.el(i);
		h = sn * g;
		g = cs * g;
		w = sqrt(h * h + f * f);
		t.el(i - 1) = w;
		cs = f / w;
		sn = h / w;
		f = x * cs + g * sn;
		g = g * cs - x * sn;
		h = y * sn;
		y = y * cs;
		for (j = 0; j < n; j++)
		    {
		    x = real(v.el(j,i - 1));
		    w = real(v.el(j,i));
		    v.el(j,i - 1) = x * cs + w * sn;
		    v.el(j,i) = w * cs - x * sn;
		    }
		w = sqrt(h * h + f * f);
		s.el(i - 1) = w;
		cs = f / w;
		sn = h / w;
		f = cs * g + sn * y;
		x = cs * y - sn * g;
		for (j = 0; j < n; j++)
		    {
		    y = real(u.el(j,i - 1));
		    w = real(u.el(j,i));
		    u.el(j,i - 1) = y * cs + w * sn;
		    u.el(j,i) = w * cs - y * sn;
		    }
		}
	    t.el(L) = 0.0;
	    t.el(k) = f;
	    s.el(k) = x;
cout << "s.el(k) is " << s.el(k) SP k << endl;
	    }
	/*
	   CONVERGENCE
	 */
	if (w >= 0.0) continue;
	s.el(k) = -w;
	for (j = 0; j < n; j++)
	    v.el(j,k) = -v.el(j,k);
	}
    /*
       SORT SINGULAR VALUES
     */
    for (k = 0; k < n; k++)	/* sort descending */
	{
	g = -1.0;
	j = k;
	for (i = k; i < n; i++)	/* sort descending */
	    {
	    if (s.el(i) <= g) continue;
	    g = s.el(i);
	    j = i;
	    }
	if (j == k) continue;
	s.el(j) = s.el(k);
	s.el(k) = g;
	for (i = 0; i < n; i++)
	    {
	    q = v.el(i,j);
	    v.el(i,j) = v.el(i,k);
	    v.el(i,k) = q;
	    }
	for (i = 0; i < n; i++)
	    {
	    q = u.el(i,j);
	    u.el(i,j) = u.el(i,k);
	    u.el(i,k) = q;
	    }
	}
    /*
       BACK TRANSFORMATION
     */
    for (k = nM1; k >= 0; k--)
	{
	if (b.el(k) == 0.0) continue;
	q = -a.el(k,k) / abs(a.el(k,k));
	for (j = 0; j < nu; j++)
	    u.el(k,j) *= q;
	for (j = 0; j < nu; j++)
	    {
	    q = 0.0;
	    for (i = k; i < m; i++)
		q += conj(a.el(i,k)) * u.el(i,j);
	    q /= abs(a.el(k,k)) * b.el(k);
	    for (i = k; i < m; i++)
		u.el(i,j) -= q * a.el(i,k);
	    }
	}
    if (n > 1)
	{
	for (k = n - 2; k >= 0; k--)
	    {
	    k1 = k + 1;
	    if (c.el(k1) == 0.0) continue;
	    q = -conj(a.el(k,k1)) / abs(a.el(k,k1));
	    for (j = 0; j < nv; j++)
		v.el(k1,j) *= q;
	    for (j = 0; j < nv; j++)
		{
		q = 0.0;
		for (i = k1; i < n; i++)
		    q += a.el(k,i) * v.el(i,j);
		q = q / (abs(a.el(k,k1)) * c.el(k1));
		for (i = k1; i < n; i++)
		    v.el(i,j) -= q * conj(a.el(k,i));
		}
	    }
	}
    } /* CSVD */

void SVDcomplex(const Matrix& Mre, const Matrix& Mim, Matrix& Ure, 
	        Matrix& Uim, Vector& d, Matrix& Vre, Matrix& Vim)
    {
    int m = Mre.Nrows(), n = Mre.Ncols();
    cout << "in SVDcomplex, m,n are " << m SP n << endl;
    ComplexMatrix U,V;
    if(m < n)	// A = u d v^dag,  so A^dag = v d u^dag = (U-Output) d (V-Output)^dag
		// my notation: A = U d V, so V = v^dag, U = u
		// So here with A^dag as input, V = (U-output)^dag, U = V-Output
	{
	Matrix mret = Mre.t(), mimt = -Mim.t();
	ComplexMatrix A(mret,mimt);
	CSVD (A, U, d, V);
	Vre = U.RealMat().t();
	Vim = -U.ImMat().t();
	Ure = V.RealMat();
	Uim = V.ImMat();
	return;
	}
    ComplexMatrix A(Mre,Mim);
    CSVD (A, U, d, V);
    Ure = U.RealMat();
    Uim = U.ImMat();
    Vre = V.RealMat().t();
    Vim = -V.ImMat().t();
    }


