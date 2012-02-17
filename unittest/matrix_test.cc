#include "test.h"
#include "matrix.h"
#include <boost/test/unit_test.hpp>
#include "boost/format.hpp"

using namespace std;
using boost::format;

struct MatrixDefaults
{
    MatrixDefaults() {} 
};

BOOST_FIXTURE_TEST_SUITE(MatrixTest,MatrixDefaults)

TEST(EigenValues)
    {
    Matrix M(2,2);
    M(1,1) = 5;
    CHECK_EQUAL(M(1,1),5);
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
    cout << format("Avg err is %.2E") % sqrt(sumerrsq/(n*m)) << endl;
    CHECK(sumerrsq < 1E-10);
    }


BOOST_AUTO_TEST_SUITE_END()

