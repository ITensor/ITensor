#include "test.h"
#include "itensor/qn.h"

using namespace itensor;
using namespace std;


TEST_CASE("QNTest")
{


SECTION("QNVal")
    {
    SECTION("Constructor")
        {
        auto q1 = QNVal(1);
        auto q2 = QNVal(2);
        auto q3 = QNVal(3);
        CHECK(q1.val() == 1);
        CHECK(q1.mod() == 1);
        CHECK(q2.val() == 2);
        CHECK(q2.mod() == 1);
        CHECK(q3.val() == 3);
        CHECK(q3.mod() == 1);

        auto qm1 = QNVal(0,2);
        CHECK(qm1.val() == 0);
        CHECK(qm1.mod() == 2);

        auto qm2 = QNVal(1,2);
        CHECK(qm2.val() == 1);
        CHECK(qm2.mod() == 2);

        auto qm3 = QNVal(2,2);
        CHECK(qm3.val() == 0);
        CHECK(qm3.mod() == 2);
        }

    SECTION("Negation")
        {
        auto q1 = QNVal(5);
        CHECK(q1.val()==5);
        auto mq1 = -q1;
        CHECK(mq1.val()==-5);

        auto q2 = QNVal(1,2);
        CHECK(q2.val()==1);
        auto mq2 = -q2;
        CHECK(mq2.val()==1);

        auto q3 = QNVal(0,2);
        CHECK(q3.val()==0);
        auto mq3 = -q3;
        CHECK(mq3.val()==0);
        }
    }

SECTION("Basic QN Constructors")
    {
    SECTION("Single Integer")
        {
        auto q0 = QN(0);
        CHECK(q0(1)==0);
        CHECK(q0.mod(1)==1);

        auto q1 = QN(0);
        CHECK(q1(1)==0);
        CHECK(q1.mod(1)==1);
        }

    SECTION("Single, Mod Factor")
        {
        auto q1 = QN({1,2});
        CHECK(q1(1)==1);
        CHECK(q1.mod(1)==2);

        auto q2 = QN({2,4});
        CHECK(q2(1)==2);
        CHECK(q2.mod(1)==4);

        auto q3 = QN({0,-1});
        CHECK(q3(1)==0);
        CHECK(q3.mod(1)==-1);
        }

    SECTION("Multiple QNVal") 
        {
        auto q1 = QN({1,2},{2,-1});
        CHECK(q1(1)==1);
        CHECK(q1.mod(1)==2);
        CHECK(q1(2)==2);
        CHECK(q1.mod(2)==-1);

        auto q2 = QN({2,1},{3,2},{0,-1});
        CHECK(q2(1)==2);
        CHECK(q2.mod(1)==1);
        CHECK(q2(2)==1);
        CHECK(q2.mod(2)==2);
        CHECK(q2(3)==0);
        CHECK(q2.mod(3)==-1);
        }

    SECTION("Multiple Integer") 
        {
        auto q1 = QN(1,2);
        CHECK(q1(1)==1);
        CHECK(q1.mod(1)==1);
        CHECK(q1(2)==2);
        CHECK(q1.mod(2)==1);

        auto q2 = QN(3,4,5);
        CHECK(q2(1)==3);
        CHECK(q2.mod(1)==1);
        CHECK(q2(2)==4);
        CHECK(q2.mod(2)==1);
        CHECK(q2(3)==5);
        CHECK(q2.mod(3)==1);
        }
    }

SECTION("Spin")
    {
    auto q0 = QN("Sz=",0);
    CHECK(q0(1) == 0);
    CHECK(Sz(q0) == 0);

    auto q1 = QN("Sz=",1);
    CHECK(q1(1) == 1);
    CHECK(Sz(q1) == 1);

    auto q2 = QN("Sz=",2);
    CHECK(q2(1) == 2);
    CHECK(Sz(q2) == 2);
    }

SECTION("Boson")
    {
    auto q0 = QN("Nb=",0);
    CHECK(q0(1) == 0);
    CHECK(Nb(q0) == 0);

    auto q1 = QN("Nb=",1);
    CHECK(q1(1) == 1);
    CHECK(Nb(q1) == 1);

    auto q2 = QN("Nb=",2);
    CHECK(q2(1) == 2);
    CHECK(Nb(q2) == 2);
    }

SECTION("SpinBoson")
    {
    auto q = QN("Sz=",0,"Nb=",0);
    CHECK(q(1) == 0);
    CHECK(q(2) == 0);
    CHECK(Sz(q) == 0);
    CHECK(Nb(q) == 0);

    q = QN("Sz=",1,"Nb=",1);
    CHECK(q(1) == 1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == 1);
    CHECK(Nb(q) == 1);

    q = QN("Sz=",-1,"Nb=",1);
    CHECK(q(1) == -1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == -1);
    CHECK(Nb(q) == 1);

    q = QN("Nb=",2,"Sz=",0);
    CHECK(q(1) == 0);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == 0);
    CHECK(Nb(q) == 2);
    }

SECTION("Fermion")
    {
    auto q = QN("Nf=",0);
    CHECK(q(1) == 0);
    CHECK(Nf(q) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q,1));

    q = QN("Nf=",1);
    CHECK(q(1) == 1);
    CHECK(Nf(q) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q,1));

    q = QN("Nf=",2);
    CHECK(q(1) == 2);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q,1));

    q = QN("Nf=",3);
    CHECK(q(1) == 3);
    CHECK(Nf(q) == 3);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q,1));
    }

SECTION("FParity")
    {
    auto q = QN("Pf=",0);
    CHECK(q(1) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q));

    q = QN("Pf=",1);
    CHECK(q(1) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q));

    q = QN("Pf=",2);
    CHECK(q(1) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q));

    q = QN("Pf=",3);
    CHECK(q(1) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q));

    q = QN("Pf=",4);
    CHECK(q(1) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q));
    }

SECTION("Electron")
    {
    auto q = QN("Sz=",0,"Nf=",0);
    CHECK(q(1) == 0);
    CHECK(q(2) == 0);
    CHECK(Sz(q) == 0);
    CHECK(Nf(q) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));
    CHECK(!isFermionic(q,1));
    CHECK(isFermionic(q,2));

    q = QN("Sz=",1,"Nf=",1);
    CHECK(q(1) == 1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == 1);
    CHECK(Nf(q) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(paritySign(q)==-1);
    CHECK(isFermionic(q));

    q = QN("Sz=",-1,"Nf=",1);
    CHECK(q(1) == -1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == -1);
    CHECK(Nf(q) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(paritySign(q)==-1);
    CHECK(isFermionic(q));

    q = QN("Sz=",0,"Nf=",2);
    CHECK(q(1) == 0);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == 0);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));

    q = QN("Sz=",2,"Nf=",2);
    CHECK(q(1) == 2);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == 2);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));

    q = QN("Sz=",-2,"Nf=",2);
    CHECK(q(1) == -2);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == -2);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));

    auto Q = QN("Sz=",+1,"Nf=",1)+QN("Sz=",-1,"Nf=",1);
    CHECK(Q == QN("Sz=",0,"Nf=",2));
    CHECK(isFermionic(Q));

    Q = QN("Sz=",0,"Nf=",2)+QN("Sz=",-1,"Nf=",1);
    CHECK(Q == QN("Sz=",-1,"Nf=",3));
    CHECK(isFermionic(Q));
    }

SECTION("Z3 Clock")
    {
    auto q = QN({0,3});
    CHECK(q(1) == 0);

    q = QN({1,3});
    CHECK(q(1) == 1);

    q = QN({2,3});
    CHECK(q(1) == 2);

    q = QN({3,3});
    CHECK(q(1) == 0);

    q = QN({-1,3});
    CHECK(q(1) == 2);

    auto Q = QN({0,3})+QN({1,3});
    CHECK(Q == QN({1,3}));

    Q = QN({0,3})+QN({2,3});
    CHECK(Q == QN({2,3}));

    Q = QN({1,3})+QN({1,3});
    CHECK(Q == QN({2,3}));

    Q = QN({1,3})+QN({2,3});
    CHECK(Q == QN({0,3}));

    Q = QN({2,3})+QN({2,3});
    CHECK(Q == QN({1,3}));
    }

}
