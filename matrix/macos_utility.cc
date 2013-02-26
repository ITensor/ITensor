// utility.cc -- Various matrix utility routines

#include "matrix.h"
#include "tarray1.h"
#include "minmax.h"
#include <math.h>
#include <fstream>

#include <Accelerate/Accelerate.h>

using namespace std;

void 
Orthog(const MatrixRef& M,int num,int numpass)	// Orthonormalize a Matrix M to num cols
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
    int lwork = max(1,4*max(n,m));
    Vector Work(lwork); Work = 0;

    Q = M.t();

    int info = 0;

    //Call lapack routine

    dgeqrf_(&m, &n, Q.Store(), &m, Tau.Store(), Work.Store(), &lwork, &info);

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
    dorgqr_(&m, &n, &tlen, Q.Store(), &m, Tau.Store(), Work.Store(), &lwork, &info);
    if(info != 0) error("Error in call to dorgqr_.");

    Q = Q.t();

    } //void QRDecomp

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
    
void 
tred2(const MatrixRef& A, Vector& D, Vector& E, Matrix& Z)
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

void 
dotranspose(Real *z, int n)
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

void 
tql2(Vector& D, Vector& E, Matrix& Z)
    {
    Real eps = FloatingPointPrecision::Epsilon ();
    int n = D.Length();
    Real *z = Z.Store();
    dotranspose(z,n);
    int l;
    for (l = 1; l < n; l++)
        E.el (l - 1) = E.el (l);
    Real b = 0.0,
         f = 0.0;
    E.el (n - 1) = 0.0;
    for (l = 0; l < n; l++)
        {
        int i, j;
        Real & dl = D.el (l);
        Real & el = E.el (l);
        Real h = eps * (fabs (dl) + fabs (el));
        if (b < h) b = h;
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

            void rotate22(Real *, Real *,Real,Real,int);
            Real *zki = z + i*n;
            Real *zki1 = zki + n;
            rotate22(zki,zki1,c,s,n);
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


void 
BackupEigenValues(const MatrixRef& A, Vector& D, Matrix& Z)
    {
    if (A.Ncols() != A.Nrows() || A.Nrows() < 1)
      _merror("BackupEigenValues: Input Matrix must be square");
    Vector E; 
    tred2(A, D, E, Z); 
    tql2(D, E, Z);
    }

//Needed for Inverse and Solve
void 
ludcmp(Matrix& a,int* indx,Real* d)
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

// Needed for Inverse and Solve
// NR, 1st edition routine lubksb
void 
lubksb(Matrix& a,int* indx,Vector& b)
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

void 
Determinant(const MatrixRef& M, Real& logdet, Real& sign)
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

// Set a Matrix of eigenvectors to the unit vectors plus a small random part
// Can be used to get an initial set of vectors for david().

void 
resetev(MatrixRef& evecs)
    { evecs.Randomize(); evecs *= 0.1; evecs += 1.0; }

// Use heapsort to sort a Vector in place. Adaptation of Num. Rec. 1st edition
// routine sort.c.

void 
Sort(Vector& V)
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


void 
EigenValues(const MatrixRef& A, Vector& D, Matrix& Z)
    {
    __CLPK_integer N = A.Ncols();
    if (N != A.Nrows() || A.Nrows() < 1)
      _merror("EigenValues: Input Matrix must be square");

    char jobz = 'V';
    char uplo = 'U';
    __CLPK_integer lwork = max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    __CLPK_doublereal work[lwork];
    __CLPK_integer info;
    
    D.ReDimension(N);
    Z = A;

    dsyev_(&jobz,&uplo,&N,Z.Store(),&N,D.Store(),work,&lwork,&info);

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

void 
GeneralizedEV(const MatrixRef& A, const MatrixRef& B, Vector& D, Matrix& Z)
    {
    __CLPK_integer N = A.Ncols();
    if (N != A.Nrows() || A.Nrows() < 1)
      _merror("EigenValues: Input Matrix must be square");

    int itype = 1; //A x = lambda B x type problem
    char jobz = 'V';
    char uplo = 'U';
    __CLPK_integer lwork = max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    __CLPK_doublereal work[lwork];
    __CLPK_integer info;
    
    D.ReDimension(N);
    Z = A;
    Matrix BB(B);//Need to copy since BB gets overwritten

    dsygv_(&itype,&jobz,&uplo,&N,Z.Store(),&N,BB.Store(),&N,D.Store(),work,&lwork,&info);

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

void GenEigenValues(const MatrixRef& A, Vector& Re, Vector& Im)
    {
    if (A.Ncols() != A.Nrows() || A.Nrows() < 1)
        _merror("GenEigenValues: Input Matrix must be square");

    if(A.Ncols() == 1)
        {
        Re = A.Column(1);
        Im = Vector(1); Im(1) = 0.0;
        return;
        }

    int n = A.Ncols();

    char jobvl = 'N'; //don't compute left evecs
    char jobvr = 'N'; //don't compute right evecs

    int lwork = max(1,4*n);
    Vector Work(lwork); Work = 0;

    Matrix Z = A;
    Re = A.Column(1);
    Im = Re;

    Vector noevecs(2); noevecs = 0; //Shouldn't be referenced by dgeev
    int num_evecs = 1;

    int info = 0;

    //Call routine
    //      JOBVL JOBVR  N   A       LDA   WR         WI         VL             LDVL       VR               LDVR      WORK        LWORK   INFO 
    dgeev_(&jobvl,&jobvr,&n,Z.Store(),&n,Re.Store(),Im.Store(),noevecs.Store(),&num_evecs,noevecs.Store(),&num_evecs,Work.Store(),&lwork,&info);
	if(info != 0)
        {
        cerr << "info is " << info << endl;
        _merror("Error condition in dsyev_.");
        }
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


void HermitianEigenvalues(const Matrix& re, const Matrix& im, Vector& evals,
	                                Matrix& revecs, Matrix& ievecs)
    {
    __CLPK_integer N = re.Ncols();
    if (N != re.Nrows() || re.Nrows() < 1)
      _merror("HermitianEigenValues: Input Matrix must be square");
    if(im.Ncols() != N || im.Nrows() != N)
      _merror("HermitianEigenValues: im not same dimensions as re");

    Matrix imt(im.t());
    ComplexMatrix H(re,imt), evecs(re,imt);

    char jobz = 'V';
    char uplo = 'U';
    __CLPK_integer lwork = max(1,3*N-1);//max(1, 1+6*N+2*N*N);
    __CLPK_doublecomplex work[lwork];
    __CLPK_doublereal rwork[lwork];
    __CLPK_integer info;
    
    evals.ReDimension(N);

    zheev_(&jobz,&uplo,&N,(__CLPK_doublecomplex*)&(H.dat[0]),&N,evals.Store(),work,&lwork,rwork,&info);
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

