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

/*
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
    CHECK(sumerrsq < 1E-10);
    }
*/

TEST(BadSVD)
    {
    int n = 4, m = 4;

    Matrix A(n,m);

    A = 0;
    A(1,3) = 3.3443453;
    A(3,3) = 3.3443453;

    //Matrix U,V; Vector D;
    //SVD(A,U,D,V);
    //cout << "D = " << endl;
    //cout << D;


    std::ifstream s("Vt");

    int nr;
    int nc;
    s.read((char*)&nr,sizeof(nr));
    s.read((char*)&nc,sizeof(nc));

    Matrix Vt(nr,nc);

    Real val;
    for(int j = 1; j <= nr; ++j)
    for(int k = 1; k <= nc; ++k)
        {
        s.read((char*)&val,sizeof(val));
        Vt(j,k) = val;
        }
    s.close();

    //cout << Vt;

    Orthog(Vt,4,2);

    //cout << Vt;
    //cout << Vt*Vt.t();

    //cout << "SVD went ok" << endl;

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

    cout << Norm(ReDiff.TreatAsVector()) << endl;
    CHECK(Norm(ReDiff.TreatAsVector()) < 1E-10);
    CHECK(Norm(ImDiff.TreatAsVector()) < 1E-10);
    }

BOOST_AUTO_TEST_SUITE_END()

