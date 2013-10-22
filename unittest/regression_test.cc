#include "test.h"
#include <boost/test/unit_test.hpp>
#include "svdalgs.h"

using namespace std;
using boost::format;

//
//The tests in this suite are
//regression tests, meaning they
//are designed to check for a specific
//bug that came up in the use of the library.
//Keeping these tests around is supposed
//to prevent the same bugs from coming up again.
//

struct RegressionDefaults
    {
    RegressionDefaults()
        {
        }

    };

BOOST_FIXTURE_TEST_SUITE(RegressionTest,RegressionDefaults)

TEST(CombinerOrder)
    {
    Index a("a",2),c("c",2);

    ITensor U(a,c);
    U(a(1),c(2)) = 1;
    //PrintDat(U);

    Combiner C(c,a);
    ITensor CU = C * U;
    //PrintDat(CU);
    ITensor UU = C * CU;
    //PrintDat(UU);

    CHECK((U-UU).norm() < 1E-10);

    Combiner D(a,c);
    CU = D * U;
    //PrintDat(CU);
    UU = D * CU;

    CHECK((U-UU).norm() < 1E-10);
    }

TEST(SVDIndexOrder)
    {
    Index a("a",2),b("b",1),c("c",2);

    ITensor z(c,a,b);
    z(c(1),a(2),b(1)) = 1;
    z(c(1),a(1),b(1)) = 1;


    ITensor V(b);

    //ITensor U(c,a); //This order succeeds
    ITensor U(a,c); //This order was failing

    ITSparse D;
    //Globals::debug2() = true;
    svd(z,U,D,V);
    //Globals::debug2() = false;

    //PrintDat(U);

    //Combiner c1(a,c),c2(c,a);

    //U = c2*U;
    //PrintDat(U);
    //U = c2*U;
    //PrintDat(U);

    ITensor nz = U*D*V;
    CHECK((z-nz).norm() < 1E-10);

    }
    //z(primed(a)(1),a(1),b(1),h(1)) = 11.8196;
    //z(primed(a)(2),a(2),b(1),h(1)) = 10.4226;
    //z(primed(a)(3),a(3),b(1),h(1)) = 9.02554;
    //z(primed(a)(1),a(1),b(2),h(1)) = 3.7886;
    //z(primed(a)(2),a(2),b(2),h(1)) = 3.11396;
    //z(primed(a)(3),a(3),b(2),h(1)) = 2.43931;
    //z(primed(a)(1),a(1),b(1),h(2)) = -1.55971;
    //z(primed(a)(3),a(3),b(1),h(2)) = 1.55971;
    //z(primed(a)(1),a(1),b(2),h(2)) = -0.753196;
    //z(primed(a)(3),a(3),b(2),h(2)) = 0.753196;
    //z(primed(a)(2),a(1),b(1),h(3)) = -1.10288;
    //z(primed(a)(3),a(2),b(1),h(3)) = -1.10288;
    //z(primed(a)(2),a(1),b(2),h(3)) = -0.53259;
    //z(primed(a)(3),a(2),b(2),h(3)) = -0.53259;
    //z(primed(a)(1),a(2),b(1),h(4)) = -1.10288;
    //z(primed(a)(1),a(2),b(2),h(4)) = -0.53259;
    //z(primed(a)(2),a(3),b(2),h(4)) = -0.53259;
    //z(primed(a)(2),a(2),b(1),h(5)) = -1.55971;
    //z(primed(a)(3),a(3),b(1),h(5)) = -1.55971;
    //z(primed(a)(1),a(1),b(2),h(5)) = -0.753196;
    //z(primed(a)(2),a(2),b(2),h(5)) = -0.753196;
    //z(primed(a)(3),a(3),b(2),h(5)) = -0.753196;

TEST(SVDArrows)
    {
    Index l("l",2),r("r",2);
    IQIndex L("L",l,QN(1,1),In),R("R",r,QN(1,1),Out);

    IQTensor AA(L,R);

    ITensor block(l,r);
    block.randomize();
    AA += block;

    const QN Zero;
    CHECK_EQUAL(div(AA),Zero);

    IQTensor U(L),V(R);
    IQTSparse D;
    svd(AA,U,D,V);

    CHECK_EQUAL(div(U),Zero);
    CHECK_EQUAL(div(V),Zero);
    }

TEST(ExpandIndex)
    {
    //
    //ITensor::expandIndex was
    //failing when the first Index 
    //had m == 1 and the second or
    //third Index was to be expanded.
    //The reason is that expandIndex
    //assumed a stable Index order
    //whereas the ITensor constructor
    //moves m==1's to the back
    //
    Index l("l");

    Index emp("emp"),occ("occ");
    IQIndex S("S",emp,QN(0,0),
                  occ,QN(1,0),
                  Out);

    ITensor oo(l,occ,primed(occ));
    oo(l(1),occ(1),primed(occ)(1)) = 1;

    oo.expandIndex(occ,S,offset(S,occ));

    oo.expandIndex(primed(occ),primed(S),offset(S,occ));

    CHECK_CLOSE(0,oo(S(1),primed(S)(1),l(1)),1E-5);
    CHECK_CLOSE(1,oo(S(2),primed(S)(2),l(1)),1E-5);
    }

TEST(ConvertToITensor)
    {
    IQIndex L("L",Index("l"),QN(),Out);
    Index emp("emp"),occ("occ");
    IQIndex S("S",emp,QN(0,0),
                  occ,QN(1,0),
                  Out);

    IQTensor T(L,S,primed(S));
    T(L(1),S(1),primed(S)(1)) = 1;
    T(L(1),S(2),primed(S)(2)) = 1;

    ITensor t = T.toITensor();

    CHECK_CLOSE(1,t(L(1),S(1),primed(S)(1)),1E-5);
    CHECK_CLOSE(1,t(L(1),S(2),primed(S)(2)),1E-5);
    CHECK_CLOSE(0,t(L(1),S(2),primed(S)(1)),1E-5);
    CHECK_CLOSE(0,t(L(1),S(2),primed(S)(1)),1E-5);
    }


/*
TEST(IndexOrder)
    {
    //
    //The ITensor contracting product code,
    //specifically the toMatrixProd part,
    //was mis-identifying certain products 
    //as matrix products when they were not.
    //

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
    */

TEST(ComplexAddition)
    {
    //EMS Oct 21, 2013
    //Bug was happening because a below has different
    //Index order from b but complex addition code 
    //did not account for this!
    Index L1("L1",6),
          S1("S1",2);

    ITensor a(L1,S1);
    a(L1(1),S1(2)) = 21;
    a(L1(2),S1(2)) = 22;
    a.scaleTo(-1.2243);

    ITensor ca = a*Complex_i;

    ITensor b(S1,L1);
    b(S1(1),L1(3)) = 13;
    b(S1(1),L1(4)) = 14;
    b.scaleTo(9.3435);

    ITensor r = b + ca;

    CHECK((realPart(r)-b).norm() < 1E-12);
    CHECK((imagPart(r)-a).norm() < 1E-12);
    }

BOOST_AUTO_TEST_SUITE_END()
