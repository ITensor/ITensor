#include "test.h"
#include "matrix.h"
#include <boost/test/unit_test.hpp>
#include "boost/format.hpp"
#include "math.h"

#include <fstream>

using namespace std;
using boost::format;

struct MatrixDefaults
    {
    int N; 
    MatrixDefaults() 
        :
        N(20)
        {} 
    };

BOOST_FIXTURE_TEST_SUITE(MatrixTest,MatrixDefaults)

TEST(MatrixVectorMultiply)
    {
    const int N = 10;

    Matrix M(N,N);
    M.Randomize();
    M *= -0.23451;

    Vector V(10);
    V.Randomize();
    V *= 1.83235;

    Vector R1 = M * V;

    const Real fac = -2.24;
    Vector R2 = (fac*M.t()) * V;

    for(int r = 1; r <= M.Nrows(); ++r)
        {
        Real val1 = 0,
             val2 = 0;
        for(int c = 1; c <= M.Ncols(); ++c)
            {
            val1 += M(r,c)*V(c);
            val2 += fac*M(c,r)*V(c);
            }
        CHECK_CLOSE(val1,R1(r),1E-5);
        CHECK_CLOSE(val2,R2(r),1E-5);
        }

    }

TEST(TestEigenValues)
    {
    Matrix A(N,N);
    A.Randomize();
    //A must be symmetric
    A += A.t();

    Matrix U;
    Vector D;
    EigenValues(A,D,U);

    for(int j = 1; j <= N; ++j)
        {
        Vector diff = D(j)*U.Column(j);
        diff -= A*U.Column(j);
        CHECK(Norm(diff) < 1E-10);
        }
    }

/*
TEST(GeneralizedEigenValues)
    {
    Matrix A(N,N);
    A.Randomize();
    A += A.t();

    vector<Vector> b(N);
    for(int j = 0; j < N; ++j)
        {
        b[j] = Vector(N);
        b[j].Randomize();
        b[j] *= 1./Norm(b[j]);
        }

    Matrix B(N,N);
    for(int i = 0; i < N; ++i)
    for(int j = 0; j < N; ++j)
        {
        B(i+1,j+1) = b[i]*b[j];
        }

    Matrix U;
    Vector D;
    GeneralizedEV(A,B,D,U);

    for(int j = 1; j <= N; ++j)
        {
        Vector diff = D(j)*B*U.Column(j);
        diff -= A*U.Column(j);
        if(Norm(diff) > 1E-8)
            {
            cerr << format("j = %d: Norm(diff) = %.3E\n") % j % Norm(diff);
            cerr << format("j = %d: D(j) = %.3E\n") % j % D(j);
            }
        CHECK(Norm(diff) < 1E-7);
        }
    }
    */

TEST(TestSVD)
    {
    int n = 200, m = 400;
    Matrix A(n,m);
    Matrix dd(n,n); dd = 0.0;
    for(int i = 1; i <= n; i++)
        {
        Real eig = pow(0.5,5*(i-1));
        dd(i,i) = eig;
        }
    Matrix uu(n,n), vv(n,m);
    uu.Randomize(); vv.Randomize();
    A = uu * dd * vv;

    Matrix U,V;  Vector D;
    SVD(A,U,D,V);

    Matrix DD(n,n); DD = 0.0; DD.Diagonal() = D;
    Matrix err = A - U * DD * V;
    Real sumerrsq = Trace(err * err.t());
    //cout << format("Avg err is %.2E") % sqrt(sumerrsq/(n*m)) << endl;
    CHECK(sumerrsq < 1E-12);
    }

TEST(TestSVDComplex)
    {
    const int n = 10,
              m = 20;
    Matrix Are(n,m),
           Aim(n,m);

    Are.Randomize();
    Aim.Randomize();

    Matrix Ure,Uim,Vre,Vim;
    Vector D;
    SVDComplex(Are,Aim,Ure,Uim,D,Vre,Vim);

    Matrix DD(D.Length(),D.Length()); DD = 0;
    for(int i = 1; i <= D.Length(); ++i) DD(i,i) = D(i);

    Matrix ReDiff = Are-(Ure*DD*Vre-Uim*DD*Vim);
    Matrix ImDiff = Aim-(Ure*DD*Vim+Uim*DD*Vre);

    CHECK(Norm(ReDiff.TreatAsVector()) < 1E-10);
    CHECK(Norm(ImDiff.TreatAsVector()) < 1E-10);

    //
    // Check nrows > ncols case
    //
    Matrix Bre(m,n),
           Bim(m,n);
    Bre.Randomize();
    Bim.Randomize();
    SVDComplex(Bre,Bim,Ure,Uim,D,Vre,Vim);
    DD = 0;
    for(int i = 1; i <= D.Length(); ++i) DD(i,i) = D(i);

    ReDiff = Bre-(Ure*DD*Vre-Uim*DD*Vim);
    ImDiff = Bim-(Ure*DD*Vim+Uim*DD*Vre);

    CHECK(Norm(ReDiff.TreatAsVector()) < 1E-10);
    CHECK(Norm(ImDiff.TreatAsVector()) < 1E-10);
    }

TEST(TestRecursiveComplexSVD)
    {
    const int n = 100,
              m = 300;
    Matrix Are(n,m),
           Aim(n,m);

    Are.Randomize();
    Aim.Randomize();

    Matrix Ure,Uim,Vre,Vim;
    Vector D;
    SVD(Are,Aim,Ure,Uim,D,Vre,Vim,1E-4);

    Matrix DD(n,n); 
    DD = 0.0; 
    DD.Diagonal() = D;

    Matrix err_re = Are - (Ure*DD*Vre - Uim*DD*Vim);
    Matrix err_im = Aim - (Ure*DD*Vim + Uim*DD*Vre);
    Real sumerrsq = Trace(err_re * err_re.t() + err_im * err_im.t());
    //cout << format("Avg err is %.2E") % sqrt(sumerrsq/(n*m)) << endl;
    CHECK(sumerrsq < 1E-12);
    }

TEST(TestHermitianEigs)
    {
    const int n = 40;
    Matrix Are(n,n),
           Aim(n,n);

    Are.Randomize();
    Aim.Randomize();
    Are = Are + Are.t();
    Aim = Aim - Aim.t();

    Matrix Ure,Uim;
    Vector D;
    HermitianEigenvalues(Are,Aim,D,Ure,Uim);

    Matrix DD(D.Length(),D.Length());
    DD = 0;
    for(int i = 1; i <= D.Length(); ++i) 
        DD(i,i) = D(i);

    Matrix ReDiff = Are-(Ure*DD*Ure.t()+Uim*DD*Uim.t());
    Matrix ImDiff = Aim-(-Ure*DD*Uim.t()+Uim*DD*Ure.t());

    CHECK(Norm(ReDiff.TreatAsVector()) < 1E-10);
    CHECK(Norm(ImDiff.TreatAsVector()) < 1E-10);
    }

TEST(TestRealDiag)
    {
    const int N = 100;
    Matrix A(N,N);

    //A(1,1) = 0.979413;
    //A(1,2) = 0.2018691;
    //A(2,1) = -0.921685;
    //A(2,2) = 0.387939;
    //A *= 1./sqrt(2);
    A.Randomize();

    //cout << "A = \n" << A << endl;

    Matrix Ure,Uim;
    Vector Dre,Dim;
    GenEigenValues(A,Dre,Dim,Ure,Uim);

    //cout << "Dre = \n" << Dre << endl;
    //cout << "Dim = \n" << Dim << endl;

    //cout << "Ure = \n" << Ure << endl;
    //cout << "Uim = \n" << Uim << endl;

    Matrix DDr(N,N),
           DDi(N,N);
    DDr = 0;
    DDi = 0;
    for(int i = 1; i <= N; ++i) 
        {
        DDr(i,i) = Dre(i);
        DDi(i,i) = Dim(i);
        }

    //Act A onto U and compare to D*U,
    //separating real and imaginary pieces

    //cout << "Re[A*U] = \n" <<  (A*Ure) << endl;
    //cout << "Re[U*D] = \n" << (Ure*DDr-Uim*DDi) << endl;

    Matrix ReDiff = A*Ure - (Ure*DDr-Uim*DDi);
    Matrix ImDiff = A*Uim - (Ure*DDi+Uim*DDr);

    //cout << (Norm(ReDiff.TreatAsVector())) << endl;
    CHECK(Norm(ReDiff.TreatAsVector()) < 1E-12);
    CHECK(Norm(ImDiff.TreatAsVector()) < 1E-12);
    }

TEST(TestNormalMatrixDiag)
    {
    //
    // For a normal matrix A such that A.t()*A == A*A.t()
    // the eigenvectors should be orthonormal such that
    // U is unitary and
    // A = U*D*U^\dagger
    //
    const int N = 3;
    Matrix A(N,N);
    A = 0;

    //A(1,1) = 1.;
    //A(1,2) = 1.;
    //A(2,1) = -1.;
    //A(2,2) = 1.;
    //A *= 1./sqrt(2);

    //
    // Example of a matrix A that is normal i.e. A.t()*A == A*A.t()
    // but is not unitary, symmetric, or anti-symmetric
    //
    A = 0;
    A(1,1) = 1.;
    A(1,2) = 1.;
    A(2,2) = 1.;
    A(2,3) = 1.;
    A(3,1) = 1.;
    A(3,3) = 1.;

    //cout << "A = \n" << A << endl;

    Matrix Ure,Uim;
    Vector Dre,Dim;
    GenEigenValues(A,Dre,Dim,Ure,Uim);

    //cout << "Dre = \n" << Dre << endl;
    //cout << "Dim = \n" << Dim << endl;

    //cout << "Ure = \n" << Ure << endl;
    //cout << "Uim = \n" << Uim << endl;

    Matrix DDr(N,N),
           DDi(N,N);
    DDr = 0;
    DDi = 0;
    for(int i = 1; i <= N; ++i) 
        {
        DDr(i,i) = Dre(i);
        DDi(i,i) = Dim(i);
        }

    Matrix ReDiff = A - (Ure*DDr*Ure.t()+Ure*DDi*Uim.t()+Uim*DDr*Uim.t()-Uim*DDi*Ure.t());
    Matrix ImPart = -Ure*DDr*Uim.t()+Ure*DDi*Ure.t()+Uim*DDr*Ure.t()+Uim*DDi*Uim.t();

    //cout << (Norm(ReDiff.TreatAsVector())) << endl;
    //cout << (Norm(ImPart.TreatAsVector())) << endl;
    CHECK(Norm(ReDiff.TreatAsVector()) < 1E-12);
    CHECK(Norm(ImPart.TreatAsVector()) < 1E-12);
    }

TEST(ComplexOrthog)
    {
    //
    // For a normal matrix A such that A.t()*A == A*A.t()
    // the eigenvectors should be orthonormal such that
    // U is unitary and
    // A = U*D*U^\dagger
    //
    const int N = 40;
    Matrix R(N,N),
           I(N,N);

    R.Randomize();
    I.Randomize();

    //cout << "R = \n" << R << endl;
    //cout << "I = \n" << I << endl;

    Orthog(R,I);

    //cout << "R = \n" << R << endl;
    //cout << "I = \n" << I << endl;

    Matrix Ore = R.t()*R + I.t()*I;
    //cout << "Ore = \n" << Ore << endl;

    Matrix Oim = R.t()*I - I.t()*R;
    //cout << "Oim = \n" << Oim << endl;

    Real re_err = 0;
    for(int r = 1; r <= N; ++r)
    for(int c = r+1; c <= N; ++c)
        {
        re_err += fabs(Ore(r,c));
        }

    for(int j = 1; j <= N; ++j)
        CHECK(fabs(Ore(j,j)-1) < 1E-12);


    CHECK(re_err < 1E-12);
    CHECK(Norm(Oim.TreatAsVector()) < 1E-12);
    }

TEST(RectQR)
    {
    const
    int r = 4,
        c = 10,
        m = min(r,c);

    Matrix M(r,c);
    M.Randomize();

    //cout << "M = \n" << M << endl;

    Matrix Q,
           R;
    QRDecomp(M,Q,R);

    Matrix I(m,m);
    I = 0;
    I.Diagonal() = 1;

    //cout << "R = \n" << R << endl;
    //cout << "Q = \n" << Q << endl;
    //cout << "Q.t()*Q = \n" << Q.t()*Q << endl;
    //cout << "Q*R-M = \n" << (Q*R-M) << endl;

    //Check that diagonal elems of R are > 0
    for(int j = 1; j <= m; ++j)
        {
        CHECK(R(j,j) > 0);
        }

    CHECK(Norm(Matrix(Q*R-M).TreatAsVector()) < 1E-14);
    CHECK(Norm(Matrix(Q.t()*Q-I).TreatAsVector()) < 1E-14);
    }

BOOST_AUTO_TEST_SUITE_END()

