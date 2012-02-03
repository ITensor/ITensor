#include "test.h"
#include "iqtensor.h"
#include <boost/test/unit_test.hpp>

using namespace std;
using boost::format;

struct RegressionDefaults
    {
    RegressionDefaults()
        {
        }

    };

BOOST_FIXTURE_TEST_SUITE(RegressionTest,RegressionDefaults)


TEST(IndexOrder)
    {
    //Globals::debug1() = true;

    Index l("l",2),r("r",2),s("s",2,Site);

    ITensor HP(l,s,r);
    HP(l(1),s(2),r(1)) = 4.30425;

    ITensor phi(l,s,r);
    phi(l(1),s(2),r(1)) = -0.341723;

    ITensor phialt(l,r,s);
    phialt.assignFrom(phi);

    //PrintDat(phi);
    //PrintDat(phialt);

    CHECK_CLOSE((phi-phialt).norm(),0,1E-3);

    ITensor order2 = HP * phi;
    //Print(order2.val0());

    //cout << endl << endl;

    ITensor order2alt = HP * phialt;
    //Print(order2alt.val0());


    CHECK_CLOSE(order2.val0(),order2alt.val0(),1E-5);


    }

    /*
    A1(l(1),s(1),m(1)) = -0.932365;
    A1(l(2),s(2),m(1)) = 0.355546;
    A1(l(3),s(3),m(1)) = -0.0654365;
    A1(l(1),s(2),m(2)) = -0.888391;
    A1(l(2),s(3),m(2)) = 0.459087;
    A1(l(1),s(3),m(3)) = 1;


    A2(r(5),m(1),s2(1)) = -0.0522478;
    A2(r(2),m(2),s2(1)) = 0.21836;
    A2(r(6),m(2),s2(1)) = -0.0248947;
    A2(r(1),m(3),s2(1)) = 0.0924649;
    A2(r(4),m(3),s2(1)) = -0.0422058;
    A2(r(7),m(3),s2(1)) = 0.00784149;
    A2(r(2),m(1),s2(2)) = -0.325766;
    A2(r(6),m(1),s2(2)) = -0.0166869;
    A2(r(1),m(2),s2(2)) = 0.384654;
    A2(r(4),m(2),s2(2)) = -0.0599322;
    A2(r(7),m(2),s2(2)) = -0.00477761;
    A2(r(3),m(3),s2(2)) = -0.0282654;
    A2(r(8),m(3),s2(2)) = -0.003913;
    A2(r(1),m(1),s2(3)) = -0.819829;
    A2(r(4),m(1),s2(3)) = -0.0328797;
    A2(r(7),m(1),s2(3)) = -0.00135719;
    A2(r(3),m(2),s2(3)) = -0.0811007;
    A2(r(8),m(2),s2(3)) = 0.00136376;
    A2(r(9),m(3),s2(3)) = -0.00242845;
    */

BOOST_AUTO_TEST_SUITE_END()
