#include "test.h"
#include "matrix.h"
#include <boost/test/unit_test.hpp>
#include "boost/format.hpp"
#include "math.h"

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


BOOST_AUTO_TEST_SUITE_END()

