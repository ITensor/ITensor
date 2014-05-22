#include "test.h"

#include "global.h"
#include "math.h"

using namespace itensor;
using namespace std;

Real
sqr(Real x) { return x*x; }

TEST_CASE("MatrixTest")
{

int N = 20;

SECTION("MatrixVectorMultiply")
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
        REQUIRE(fabs(val1-R1(r)) < 1E-5);
        REQUIRE(fabs(val2-R2(r)) < 1E-5);
        }

    }

SECTION("TestEigenValues")
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
        REQUIRE(Norm(diff) < 1E-10);
        }
    }


SECTION("TestSVD")
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
    REQUIRE(sumerrsq < 1E-12);
    }

SECTION("TestSVDComplex")
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

    REQUIRE(Norm(ReDiff.TreatAsVector()) < 1E-10);
    REQUIRE(Norm(ImDiff.TreatAsVector()) < 1E-10);

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

    REQUIRE(Norm(ReDiff.TreatAsVector()) < 1E-10);
    REQUIRE(Norm(ImDiff.TreatAsVector()) < 1E-10);
    }

SECTION("TestRecursiveComplexSVD")
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
    REQUIRE(sumerrsq < 1E-12);
    }

SECTION("TestHermitianEigs")
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

    REQUIRE(Norm(ReDiff.TreatAsVector()) < 1E-10);
    REQUIRE(Norm(ImDiff.TreatAsVector()) < 1E-10);
    }

SECTION("TestRealDiag")
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
    REQUIRE(Norm(ReDiff.TreatAsVector()) < 1E-12);
    REQUIRE(Norm(ImDiff.TreatAsVector()) < 1E-12);
    }

SECTION("TestNormalMatrixDiag")
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
    REQUIRE(Norm(ReDiff.TreatAsVector()) < 1E-12);
    REQUIRE(Norm(ImPart.TreatAsVector()) < 1E-12);
    }

SECTION("ComplexOrthog")
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
        REQUIRE(fabs(Ore(j,j)-1) < 1E-12);


    REQUIRE(re_err < 1E-12);
    REQUIRE(Norm(Oim.TreatAsVector()) < 1E-12);
    }

SECTION("RectQR")
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
    //for(int j = 1; j <= m; ++j)
    //    {
    //    REQUIRE(R(j,j) > 0);
    //    if(R(j,j) <= 0) printfln("R(%d,%d) = %.5E"),j,j, R(j,j) << endl;
    //    }

    REQUIRE(Norm(Matrix(Q*R-M).TreatAsVector()) < 1E-14);
    REQUIRE(Norm(Matrix(Q.t()*Q-I).TreatAsVector()) < 1E-14);
    }

SECTION("ComplexEV")
    {
    const int N = 40;
    Matrix R(N,N),
           I(N,N);

    R.Randomize();
    I.Randomize();

    Matrix VR(N,N),
           VI(N,N);

    Vector eR(N),
           eI(N);

    ComplexEigenvalues(R,I,eR,eI,VR,VI);

    Matrix PR = R*VR - I*VI,
           PI = R*VI + I*VR;

    for(int j = 1; j <= N; ++j)
        {
        //jth eigenvector
        VectorRef vR = VR.Column(j),
                  vI = VI.Column(j);

        //Matrix times eigenvector
        Vector pR = R*vR-I*vI,
               pI = R*vI+I*vR;

        //Eigenvalue time eigenvector
        Vector zR = eR(j)*vR - eI(j)*vI,
               zI = eR(j)*vI + eI(j)*vR;

        //Should be the same
        Real nrm = 0;
        nrm += sqr(Norm(pR-zR));
        nrm += sqr(Norm(pI-zI));
        nrm = sqrt(nrm);

        REQUIRE(nrm < 1E-12);
        }
    }
}


