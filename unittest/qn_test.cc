#include "test.h"
#include "itensor/qn.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using namespace std;


TEST_CASE("QNTest")
{


SECTION("QNum")
    {
    SECTION("Constructor")
        {
        auto q1 = QNum(1);
        auto q2 = QNum(2);
        auto q3 = QNum(3);
        CHECK(q1.val() == 1);
        CHECK(q1.mod() == 1);
        CHECK(q2.val() == 2);
        CHECK(q2.mod() == 1);
        CHECK(q3.val() == 3);
        CHECK(q3.mod() == 1);

        auto qm1 = QNum(0,2);
        CHECK(qm1.val() == 0);
        CHECK(qm1.mod() == 2);

        auto qm2 = QNum(1,2);
        CHECK(qm2.val() == 1);
        CHECK(qm2.mod() == 2);

        auto qm3 = QNum(2,2);
        CHECK(qm3.val() == 0);
        CHECK(qm3.mod() == 2);
        }

    SECTION("Negation")
        {
        auto q1 = QNum(5);
        CHECK(q1.val()==5);
        auto mq1 = -q1;
        CHECK(mq1.val()==-5);

        auto q2 = QNum(1,2);
        CHECK(q2.val()==1);
        auto mq2 = -q2;
        CHECK(mq2.val()==1);

        auto q3 = QNum(0,2);
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
        CHECK(q0.val(1)==0);
        CHECK(q0.mod(1)==1);

        auto q1 = QN(0);
        CHECK(q1.val(1)==0);
        CHECK(q1.mod(1)==1);
        }

    SECTION("Single, Mod Factor")
        {
        auto q1 = QN({1,2});
        CHECK(q1.val(1)==1);
        CHECK(q1.mod(1)==2);
        CHECK(q1.name(1)=="");

        auto q2 = QN({2,4});
        CHECK(q2.val(1)==2);
        CHECK(q2.mod(1)==4);
        CHECK(q2.name(1)=="");


        auto q3 = QN({0,-1});
        CHECK(q3.val(1)==0);
        CHECK(q3.mod(1)==-1);
        CHECK(q3.name(1)=="");
        }

    SECTION("Multiple QNum") 
        {
        auto q1 = QN({"A",1,2},{"B",2,-1});
        CHECK(q1.val(1)==1);
        CHECK(q1.mod(1)==2);
        CHECK(q1.name(1)=="A");
        CHECK(q1.val(2)==2);
        CHECK(q1.mod(2)==-1);
        CHECK(q1.name(2)=="B");

        auto q2 = QN({"AA",2,1},{"BB",3,2},{"CC",0,-1});
        CHECK(q2.val(1)==2);
        CHECK(q2.mod(1)==1);
        CHECK(q2.name(1)=="AA");
        CHECK(q2.val(2)==1);
        CHECK(q2.mod(2)==2);
        CHECK(q2.name(2)=="BB");
        CHECK(q2.val(3)==0);
        CHECK(q2.mod(3)==-1);
        CHECK(q2.name(3)=="CC");
        }
    }

SECTION("Spin")
    {
    auto q0 = QN({"Sz",0});
    CHECK(q0.val("Sz") == 0);
    CHECK(q0.mod("Sz") == 1);

    auto q1 = QN({"Sz",1});
    CHECK(q1.val("Sz") == 1);

    auto qn1 = QN({"Sz",-1});
    CHECK(qn1.val("Sz") == -1);

    auto q2 = QN({"Sz",2});
    CHECK(q2.val("Sz") == 2);
    }

SECTION("Boson")
    {
    auto q0 = QN({"Nb",0});
    CHECK(q0.val("Nb") == 0);
    CHECK(q0.mod("Nb") == 1);

    auto q1 = QN({"Nb",1});
    CHECK(q1.val("Nb") == 1);
    CHECK(q1.mod("Nb") == 1);

    auto q2 = QN({"Nb",2});
    CHECK(q2.val("Nb") == 2);
    CHECK(q2.mod("Nb") == 1);
    }

SECTION("SpinBoson")
    {
    auto q = QN({"Sz",0},{"Nb",0});
    CHECK(q.val("Nb") == 0);
    CHECK(q.val("Sz") == 0);

    q = QN({"Sz",1},{"Nb",0});
    CHECK(q.val("Nb") == 0);
    CHECK(q.val("Sz") == 1);

    q = QN({"Sz",1},{"Nb",1});
    CHECK(q.val("Nb") == 1);
    CHECK(q.val("Sz") == 1);

    q = QN({"Sz",-1},{"Nb",1});
    CHECK(q.val("Nb") == 1);
    CHECK(q.val("Sz") == -1);

    q = QN({"Sz",0},{"Nb",2});
    CHECK(q.val("Nb") == 2);
    CHECK(q.val("Sz") == 0);
    }

SECTION("Fermion")
    {
    auto q = QN({"Nf",0,-1});
    CHECK(q.val("Nf") == 0);
    CHECK(q.mod("Nf") == -1);
    CHECK(isFermionic(q,"Nf"));

    q = QN({"Nf",1,-1});
    CHECK(q.val("Nf") == 1);
    CHECK(q.mod("Nf") == -1);
    CHECK(isFermionic(q,"Nf"));

    q = QN({"Nf",2,-1});
    CHECK(q.val("Nf") == 2);
    CHECK(q.mod("Nf") == -1);
    CHECK(isFermionic(q,"Nf"));

    q = QN({"Nf",3,-1});
    CHECK(q.val("Nf") == 3);
    CHECK(q.mod("Nf") == -1);
    CHECK(isFermionic(q,"Nf"));
    }

SECTION("Parity, Fermionic")
    {
    auto q = QN({"Nfp",0,-2});
    CHECK(q.val("Nfp") == 0);
    CHECK(q.mod("Nfp") == -2);
    CHECK(isFermionic(q,"Nfp"));

    q = QN({"Nfp",1,-2});
    CHECK(q.val("Nfp") == 1);
    CHECK(q.mod("Nfp") == -2);
    CHECK(isFermionic(q,"Nfp"));

    q = QN({"Nfp",2,-2});
    CHECK(q.val("Nfp") == 0);
    CHECK(q.mod("Nfp") == -2);
    CHECK(isFermionic(q,"Nfp"));
    }

SECTION("Z3 Clock")
    {
    auto q = QN({"T",0,3});
    CHECK(q.val("T") == 0);

    q = QN({"T",1,3});
    CHECK(q.val("T") == 1);

    q = QN({"T",2,3});
    CHECK(q.val("T") == 2);

    q = QN({"T",3,3});
    CHECK(q.val("T") == 0);

    q = QN({"T",4,3});
    CHECK(q.val("T") == 1);

    q = QN({"T",-1,3});
    CHECK(q.val("T") == 2);

    auto Q = QN({"T",0,3})+QN({"T",1,3});
    CHECK(Q == QN({"T",1,3}));

    Q = QN({"T",1,3})+QN({"T",1,3});
    CHECK(Q == QN({"T",2,3}));

    Q = QN({"T",2,3})+QN({"T",1,3});
    CHECK(Q == QN({"T",0,3}));

    Q = QN({"T",2,3})-QN({"T",1,3});
    CHECK(Q == QN({"T",1,3}));
    }

}
