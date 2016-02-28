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

SECTION("Spin")
    {
    auto q0 = spin(0);
    CHECK(q0(1) == 0);
    CHECK(Sz(q0) == 0);

    auto q1 = spin(1);
    CHECK(q1(1) == 1);
    CHECK(Sz(q1) == 1);

    auto q2 = spin(2);
    CHECK(q2(1) == 2);
    CHECK(Sz(q2) == 2);
    }

SECTION("Boson")
    {
    auto q0 = boson(0);
    CHECK(q0(1) == 0);
    CHECK(Nb(q0) == 0);

    auto q1 = boson(1);
    CHECK(q1(1) == 1);
    CHECK(Nb(q1) == 1);

    auto q2 = boson(2);
    CHECK(q2(1) == 2);
    CHECK(Nb(q2) == 2);
    }

SECTION("SpinBoson")
    {
    auto q = spinboson(0,0);
    CHECK(q(1) == 0);
    CHECK(q(2) == 0);
    CHECK(Sz(q) == 0);
    CHECK(Nb(q) == 0);

    q = spinboson(1,1);
    CHECK(q(1) == 1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == 1);
    CHECK(Nb(q) == 1);

    q = spinboson(-1,1);
    CHECK(q(1) == -1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == -1);
    CHECK(Nb(q) == 1);

    q = spinboson(0,2);
    CHECK(q(1) == 0);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == 0);
    CHECK(Nb(q) == 2);
    }

SECTION("Fermion")
    {
    auto q = fermion(0);
    CHECK(q(1) == 0);
    CHECK(Nf(q) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q,1));

    q = fermion(1);
    CHECK(q(1) == 1);
    CHECK(Nf(q) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q,1));

    q = fermion(2);
    CHECK(q(1) == 2);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q,1));

    q = fermion(3);
    CHECK(q(1) == 3);
    CHECK(Nf(q) == 3);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q,1));
    }

SECTION("FParity")
    {
    auto q = fparity(0);
    CHECK(q(1) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q));

    q = fparity(1);
    CHECK(q(1) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(isFermionic(q));

    q = fparity(2);
    CHECK(q(1) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(isFermionic(q));
    }

SECTION("Electron")
    {
    auto q = electron(0,0);
    CHECK(q(1) == 0);
    CHECK(q(2) == 0);
    CHECK(Sz(q) == 0);
    CHECK(Nf(q) == 0);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));
    CHECK(!isFermionic(q,1));
    CHECK(isFermionic(q,2));

    q = electron(1,1);
    CHECK(q(1) == 1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == 1);
    CHECK(Nf(q) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(paritySign(q)==-1);
    CHECK(isFermionic(q));

    q = electron(-1,1);
    CHECK(q(1) == -1);
    CHECK(q(2) == 1);
    CHECK(Sz(q) == -1);
    CHECK(Nf(q) == 1);
    CHECK(Nfp(q) == 1);
    CHECK(paritySign(q)==-1);
    CHECK(isFermionic(q));

    q = electron(0,2);
    CHECK(q(1) == 0);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == 0);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));

    q = electron(2,2);
    CHECK(q(1) == 2);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == 2);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));

    q = electron(-2,2);
    CHECK(q(1) == -2);
    CHECK(q(2) == 2);
    CHECK(Sz(q) == -2);
    CHECK(Nf(q) == 2);
    CHECK(Nfp(q) == 0);
    CHECK(paritySign(q)==+1);
    CHECK(isFermionic(q));

    auto Q = electron(+1,1)+electron(-1,1);
    CHECK(Q == electron(0,2));
    CHECK(isFermionic(Q));

    Q = electron(0,2)+electron(-1,1);
    CHECK(Q == electron(-1,3));
    CHECK(isFermionic(Q));
    }

SECTION("Z3 Clock")
    {
    auto q = clock(0,3);
    CHECK(q(1) == 0);

    q = clock(1,3);
    CHECK(q(1) == 1);

    q = clock(2,3);
    CHECK(q(1) == 2);

    q = clock(3,3);
    CHECK(q(1) == 0);

    q = clock(-1,3);
    CHECK(q(1) == 2);

    auto Q = clock(0,3)+clock(1,3);
    CHECK(Q == clock(1,3));

    Q = clock(0,3)+clock(2,3);
    CHECK(Q == clock(2,3));

    Q = clock(1,3)+clock(1,3);
    CHECK(Q == clock(2,3));

    Q = clock(1,3)+clock(2,3);
    CHECK(Q == clock(0,3));

    Q = clock(2,3)+clock(2,3);
    CHECK(Q == clock(1,3));
    }
}
