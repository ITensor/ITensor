// test.cc -- Test the matrix package
#define THIS_IS_MAIN

#include "matrix.h"
#include <math.h>
#include <fstream>
#include <iomanip>

extern Real svdtruncerr;
#ifdef THIS_IS_MAIN
Real svdtruncerr = 0.0;
#endif
void svdcut(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V, Real truncerr, int maxm)
    {
    int m = A.Nrows(), n = A.Ncols();	// want m >= n
    if(m < n)
	{
	Matrix Ut, Vt;
	svdcut(A.t(),Ut,d,Vt,truncerr,maxm);
	U = Vt.t();
	V = Ut.t();
	return;
	}
    Matrix AA(A), evecs;
    Vector evals;
    Matrix Asq = AA.t() * AA;
    Real sca = 0.0;
    for(int i = 1; i <= n; i++)
	for(int j = 1; j <= n; j++)
	    if(fabs(Asq(i,j)) > sca) sca = fabs(Asq(i,j));
    if(sca > 0.0) Asq *= -1.0/sca;
    EigenValues(Asq,evals,evecs);
    evals *= -sca;
    int k = n;
    Real err = 0.0;
    while(err+evals(k) < truncerr*evals(1) && k > 1)
	err += evals(k--);
    svdtruncerr = err/evals(1);
    cout << "k is " << k << endl;
    V = evecs.Columns(1,k).t();
    U = AA * V.t();
    d.ReDimension(k);
    for(int i = 1; i <= k; i++)
	d(i) = sqrt(evals(i)),
	U.Column(i) *= 1.0 / d(i);
    }
/*
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
    */


Real ran1();
int main()
{
    cout << "sizeof(bool) is " << sizeof(bool) << endl;
    if(0)
	{
        int nrows,ncols;
        ifstream file;
        file.open("bad_matrix");
        file.read((char*)&nrows,sizeof(nrows));
        file.read((char*)&ncols,sizeof(ncols));
        cout << "nrows, ncols are " << nrows SP ncols << endl;
        Matrix Are(nrows,ncols), Aim(nrows,ncols);
        file.read((char*)Are.Store(),sizeof(Real)*nrows*ncols); 
        file.read((char*)Aim.Store(),sizeof(Real)*nrows*ncols); 

        file.read((char*)&nrows,sizeof(nrows));
        file.read((char*)&ncols,sizeof(ncols));
        cout << "nrows, ncols are " << nrows SP ncols << endl;
        Matrix Ure(nrows,ncols), Uim(nrows,ncols);
        file.read((char*)Ure.Store(),sizeof(Real)*nrows*ncols); 
        file.read((char*)Uim.Store(),sizeof(Real)*nrows*ncols); 

        file.read((char*)&nrows,sizeof(nrows));
        file.read((char*)&ncols,sizeof(ncols));
        cout << "nrows, ncols are " << nrows SP ncols << endl;
        Matrix Vre(nrows,ncols), Vim(nrows,ncols);
        file.read((char*)Vre.Store(),sizeof(Real)*nrows*ncols); 
        file.read((char*)Vim.Store(),sizeof(Real)*nrows*ncols); 

        file.read((char*)&nrows,sizeof(nrows));
        Vector d(nrows);
        file.read((char*)d.Store(),sizeof(Real)*nrows); 
        file.close();
        cout << "d length is " << nrows << endl;
        cout << "norm of Are is " << Norm(Are.TreatAsVector()) << endl;
        cout << "norm of Aim is " << Norm(Aim.TreatAsVector()) << endl;
        cout << "norm of Ure is " << Norm(Ure.TreatAsVector()) << endl;
        cout << "norm of Uim is " << Norm(Uim.TreatAsVector()) << endl;
        cout << "norm of Vre is " << Norm(Vre.TreatAsVector()) << endl;
        cout << "norm of Vim is " << Norm(Vim.TreatAsVector()) << endl;

        Real err = check_complex_SVD( Are,  Aim,  Ure, Uim,  d,  Vre,  Vim);
        cout << "err is " << err << endl;
        cout << "d = " << setprecision(15) << d;

        newcomplexSVD(Are, Aim, Ure, Uim, d, Vre, Vim);
        cout << "d = " << setprecision(15) << d;
        //cout << "Ure = " << Ure;
        //cout << "Uim = " << Uim;
        //cout << "Vre = " << Vre;
        //cout << "Vim = " << Vim;
        err = check_complex_SVD( Are,  Aim,  Ure, Uim,  d,  Vre,  Vim);
        cout << "err is " << err << endl;
        exit(0);
	}
    if(0)
	{
        int m,n;
        cin >> m >> n;
        if(1)
            {
            Matrix Mre(m,n), Mim(m,n),Ure,Uim,Vre,Vim;
            Vector d;
            for(int i = 1; i <= m; i++)
            for(int j = 1; j <= n; j++)
                Mre(i,j) = ran1(), Mim(i,j) = ran1();
            newcomplexSVD(Mre, Mim, Ure, Uim, d, Vre, Vim);
            cout << "d = " << d;
            cout << "Ure = " << Ure;
            cout << "Uim = " << Uim;
            cout << "Vre = " << Vre;
            cout << "Vim = " << Vim;
            Real err = check_complex_SVD( Mre,  Mim,  Ure, Uim,  d,  Vre,  Vim);
            cout << "err is " << err << endl;
            exit(0);
            }
	}
    if(0)
	{
        Matrix A(50,100); A.Randomize();
        //cout << A;
        Matrix U,V; Vector d;
        //SVD(A,U,d,V);
        svdcut(A,U,d,V,1.0e-6,50);
        cout << "svdtruncerr is " << svdtruncerr << endl;
        Matrix D(d.Length(),d.Length()); D = 0;
        D.Diagonal() = d;
        //cout << A - Matrix(U * D * V);
        Matrix check(A - Matrix(U * D * V));
        cout << "error is " << Norm(check.TreatAsVector()) << endl;
        cout << d;
        exit(0);
        A = A.t();
        SVD(A,U,d,V);
        D.ReDimension(d.Length(),d.Length()); D = 0;
        D.Diagonal() = d;
        cout << A - Matrix(U * D * V);
        check = A - Matrix(U * D * V);
        cout << "error is " << Norm(check.TreatAsVector()) << endl;
        cout << d;
        exit(0);
	}
    int gig = 1 << 31;
    cout << "1 << 31 is " << gig << endl;
    cout << "-(1 << 31) is " << -gig << endl;
    cout << "-(1 << 31) - 1 is " << -gig-1 << endl;
    cout << "(1 << 31) + 1 is " << gig+1 << endl;
    gig = ((long)1 << 31) - 1;
    cout << "1 << 31 - 1 is " << gig << endl;
    cout << "-(1 << 31 - 1) is " << -gig << endl;
    cout << "-(1 << 31-1) - 1 is " << -gig-1 << endl;
    cout << "(1 << 31-1) + 1 is " << gig+1 << endl;
    cout << "sizeof(short) is " << sizeof(short) << endl;
    cout << "sizeof(int) is " << sizeof(int) << endl;
    cout << "sizeof(long) is " << sizeof(long) << endl;
    cout << "sizeof(long long) is " << sizeof(long long) << endl;
    cout << "sizeof(long double) is " << sizeof(long double) << endl;
    cout << "sizeof(void *) is " << sizeof(void *) << endl;
    cout << "sizeof(size_t) is " << sizeof(size_t) << endl;
    Matrix A(6, 6);
    A = 1;
    Matrix B = Transpose(A);
    Matrix C = A * B;
    cout << "A,B,C:\n" << A << B << C;
    int i, j;

    for (i = 1; i <= C.Nrows(); i++)
	for (j = 1; j <= C.Ncols(); j++)
	    C(i, j) = j * j * j - 8 * j * j + (i * i - 5 * i) * C.Ncols() + (i == j);

    cout << "C:\n" << (Matrix)(C * A) << (Matrix)(A * C);
    Matrix DD(C);
    DD += 3.0 * DD;
    DD -= 4.0 * C;
    cout << "zero should be " << Norm(DD.TreatAsVector()) << endl;
    Matrix Ci = Inverse(C);
    cout << "Inverse(C):\n" << Ci << "C^-1*C:\n" << Ci * C;
    cout << "dot of diags: C, Ci = " << C.Diagonal() * Ci.Diagonal() << endl;
    Real res = 0.0;
    for (i = 1; i <= C.Nrows(); i++)
	res += C(i, i) * Ci(i, i);
    cout << "should be " << res << endl;
    cout << "C.t() is \n" << C.t();
    Matrix ct = C.t();
    cout << "same is \n" << ct;
    Matrix ctx = C + ct;
    cout << "C + C is \n" << C + C;
    Matrix c2 = C;
    c2 *= 2;
    cout << "2*C is \n" << c2;
    cout << "C + C.t() is (should be symm)\n" << C + C.t();
    Matrix ctt = C;
    ctt += C.t();
    cout << "same(+=) is \n" << ctt;
//    ctt=2*C; ctt += 2.0 * C.t();
//    cout << "same(2/2) is \n" << ctt/2.0;

    Orthog(C);
    cout << "orthoged C:\n" << C;
    Matrix Ct = C.t();
    cout << "C.t(), one:\n" << Ct << Ct * C << C.t() * C;

    Matrix D = A + (B * C + Inverse(C));
    cout << "D =\n" << D;
    cout << "one:\n" << (D - A - B * C) * C;

    Matrix Cpart = C.SubMatrix(1, 3, 1, 5);
    cout << "Cpart, Cpart.t() = \n" << Cpart << Cpart.t();
    Vector vp5 = C.Column(6).SubVector(1, 5);
    Vector vp3 = C.Row(5).SubVector(2, 4);
    cout << "vp5, vp3 = \n" << vp5 << vp3;
    cout << "Cpart * vp5 = , vp3 * Cpart = \n" << Cpart * vp5 <<
	vp3 * Cpart;

    Vector VC = C.Column(1) * 5;
    Real x = -C.Column(1) * C.Column(2);
    cout << "C,VC:\n" << C << VC << "x:" << x << endl;
    cout << "C * VC = \n" << C * VC;
    cout << "VC * C = \n" << VC * C;
    C.Column(1) = 1;
    Vector V3(6);
    V3 = 0;
    V3(2) = 1;
    C.Column(1) += V3 * 2;
    B = A * C;
    A = B * 2;
    cout << "A,B,C:\n" << A << B << C;
    Orthog(A);
    cout << "orth'd A:\n" << A;
    cout << "One:\n" << A.t() * A;
    cout << "Smaller one:" << Transpose(A.Columns(1, 3)) * A.Columns(1, 3);
    Matrix Asmall = A.Columns(1, 1);
    cout << "Asmall is:\n" << Asmall;
    Matrix as2 = Asmall;
    Matrix rsmall(1, 1);
// mult(as2,Asmall,rsmall,1,0);
    rsmall = as2.t() * Asmall;
    cout << "Another small one:\n" << rsmall;
// mult(as2,Asmall,A,0,1);
    A = as2 * Asmall.t();
    cout << "Big A:\n" << A;

    A = 0;
    for (i = 1; i <= C.Nrows(); i++)
	for (j = 1; j <= C.Ncols(); j++)
	    if (i == j - 1 || i == j + 1)
		A(i, j) = -1;

    Matrix Ai = Inverse(A);
    cout << "A, Ai, Ai*A:\n" << A << Ai << Ai * A;
    cout << "Total Storage in StoreLinks is now " <<
	StoreLink::TotalStorage() << endl;

    cout << "Checking EigenValues. " << endl;
    Vector E;
    EigenValues(A, E, B);
    cout << "A, B, E's:\n" << A << B << E;
    Vector eigcheck = A * B.Column(1) - E(1) * B.Column(1);
    cout << "Check on first col eig vec: " << Norm(eigcheck) << endl;
    Matrix Bt = B.t();
    cout << ".t():\n" << Bt << "One:\n" << Bt * B;
    Orthog(B);
    Bt = Transpose(B);
    cout << "B after orthog:\n" << B << Bt;
    C = Bt * B;
    cout << "One:\n" << C;
    Matrix BS = B.SubMatrix(1,5,2,3);
    Matrix BSS = B.SubMatrix(2,4,2,3);
    cout << "BS * BSS.t() = " << BS * BSS.t();

    Matrix Bi = Inverse(B);
    cout << "B, Bi, Bi*B, B*Bi:\n" << B << Bi << Bi * B << B * Bi;
    Bi = 1;
    cout << "Bi from Solve:\n" << Solve(B, Bi);

    Vector V(6);
    V = 1;
    V *= 5;
    Vector V2;
    V2 = V;
    cout << "V, V2:\n" << V << V2;

    cout << "B, V * B: " << B << V * B;

    Matrix P(8, 4);
    for (i = 1; i <= P.Nrows(); i++)
	for (j = 1; j <= P.Ncols(); j++)
	    P(i, j) = 1;
    P += 2;
    cout << "Matrix plus a scalar:\n" << P << "plus another:\n" << 1 + P;
    VectorRefNoLink vv;
    vv << V;
    vv = V2;

/*
Vector in(128),outre(65),outim(65);
in = 0.0;
in(1) = 1.0;
FFT(in,outre,outim);
cout << "FFT test: Real out should be const:\n" << outre;
cout << "FFT test: Imag out should be 0:\n" << outim;
  
for (i = 0 ; i < 128 ; i++)
    in[i] = sin(i*M_PI/16);
FFT(in,outre,outim);
cout << "FFT test: FT of a sine wave\n";
cout << "Real out is :\n" << outre;
cout << "Imag out is :\n" << outim;
*/
    if(1)
	{
	Matrix A(50,100); A.Randomize();
	cout << A;
	Matrix U,V; Vector d;
	SVD(A,U,d,V);
	Matrix D(d.Length(),d.Length()); D = 0;
	D.Diagonal() = d;
	cout << A - Matrix(U * D * V);
	Matrix check(A - Matrix(U * D * V));
	cout << "error is " << Norm(check.TreatAsVector()) << endl;
	cout << d;
	A = A.t();
	SVD(A,U,d,V);
	D.ReDimension(d.Length(),d.Length()); D = 0;
	D.Diagonal() = d;
	cout << A - Matrix(U * D * V);
	check = A - Matrix(U * D * V);
	cout << "error is " << Norm(check.TreatAsVector()) << endl;
	cout << d;
	}
    if(1)
	{
	Matrix A(50,100); A.Randomize();
	cout << A;
	Matrix U,V; Vector d;
void newSVD(const MatrixRef& A, Matrix& U, Vector& d, Matrix& V);
	newSVD(A,U,d,V);
	Matrix D(d.Length(),d.Length()); D = 0;
	D.Diagonal() = d;
	cout << A - Matrix(U * D * V);
	Matrix check(A - Matrix(U * D * V));
	cout << "new SVD error is " << Norm(check.TreatAsVector()) << endl;
	cout << d;
	A = A.t();
	newSVD(A,U,d,V);
	D.ReDimension(d.Length(),d.Length()); D = 0;
	D.Diagonal() = d;
	cout << A - Matrix(U * D * V);
	check = A - Matrix(U * D * V);
	cout << "new SVD error is " << Norm(check.TreatAsVector()) << endl;
	cout << d;
	}
    if(1)
	{
	Matrix hre(2,2), him(2,2);
	him = hre = 0.0;
	hre(2,2) = 1.0;
	him(1,2) = -1.0;
	him(2,1) = 1.0;
	Matrix revecs,ievecs;
	Vector evals;
	HermitianEigenvalues(hre,him,evals,revecs,ievecs);
	cout << "evals are " << evals;
	cout << "revecs are " << revecs;
	cout << "ievecs are " << ievecs;
	cout << "exact from Mathematica is \n{{1.61803, -0.618034}, {{0. - 0.525731 I, \n 0.850651}, {0. - 0.850651 I, -0.525731}}}\n";
	}

    }
