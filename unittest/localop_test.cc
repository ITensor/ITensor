#include "test.h"
#include "itensor/mps/localop.h"
#include "itensor/mps/localmpo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

TEST_CASE("LocalOp")
{
auto s1 = Index(1,"Site");
auto s2 = Index(2,"Site");
auto h0 = Index(4,"Link");
auto h1 = Index(4,"Link");
auto h2 = Index(4,"Link");
auto l0 = Index(10,"Link");
auto l1 = Index(10,"Link");
auto l2 = Index(10,"Link");

auto S1 = Index(QN(-1),1,
                QN(+1),1,
                "Site");
auto S2 = Index(QN(-1),1,
                QN(+1),1,
                "Site");
auto H0 = Index(QN(-2),4,
                QN(+0),8,
                QN(+2),4,
                "Link");
auto H1 = Index(QN(-2),4,
                QN(+0),8,
                QN(+2),4,
                "Link");
auto H2 = Index(QN(-2),4,
                QN(+0),8,
                QN(+2),4,
                "Link");
auto L0 = Index(QN(-2),4,
                QN(+0),8,
                QN(+2),4,
                "Link");
auto L1 = Index(QN(-2),4,
                QN(+0),8,
                QN(+2),4,
                "Link");
auto L2 = Index(QN(-2),4,
                QN(+0),8,
                QN(+2),4,
                "Link");


SECTION("Product")
    {
    SECTION("Bulk Case - 2 center site")
        {
        auto Op1 = randomITensor(s1,prime(s1),h0,h1);
        auto Op2 = randomITensor(s2,prime(s2),h1,h2);
        auto L = randomITensor(l0,prime(l0),h0);
        auto R = randomITensor(l2,prime(l2),h2);
        auto lop = LocalOp(Op1,Op2,L,R);
        auto psi = randomITensor(l0,s1,s2,l2);
        auto Hpsi = ITensor();
        lop.product(psi,Hpsi);
        CHECK(hasIndex(Hpsi,s1));
        CHECK(hasIndex(Hpsi,s2));
        CHECK(hasIndex(Hpsi,l0));
        CHECK(hasIndex(Hpsi,l2));
        }

    SECTION("Bulk Case - 1 center site")
        {
        auto Op1 = randomITensor(s1,prime(s1),h0,h1);
        auto L = randomITensor(l0,prime(l0),h0);
        auto R = randomITensor(l1,prime(l1),h1);
        auto lop = LocalOp(Op1,L,R);
        auto psi = randomITensor(l0,s1,l1);
        auto Hpsi = ITensor();
        lop.product(psi,Hpsi);
        CHECK(hasIndex(Hpsi,s1));
        CHECK(hasIndex(Hpsi,l0));
        CHECK(hasIndex(Hpsi,l1));
        }
    
    SECTION("Bulk Case - 0 center site")
        {
        auto L = randomITensor(l0,prime(l0),h0);
        auto R = randomITensor(l1,prime(l1),h0);
        auto lop = LocalOp(L,R);
        auto psi = randomITensor(l0,l1);
        auto Hpsi = ITensor();
        lop.product(psi,Hpsi);
        CHECK(hasIndex(Hpsi,l0));
        CHECK(hasIndex(Hpsi,l1));
        }

    }

SECTION("Diag")
    {
    SECTION("Bulk Case - ITensor")
        {
        auto Op1 = randomITensor(s1,prime(s1),h0,h1);
        auto Op2 = randomITensor(s2,prime(s2),h1,h2);
        auto L = randomITensor(l0,prime(l0),h0);
        auto R = randomITensor(l2,prime(l2),h2);
        auto lop = LocalOp(Op1,Op2,L,R);
        auto diag = lop.diag();
        CHECK(hasIndex(diag,s1));
        CHECK(hasIndex(diag,s2));
        CHECK(hasIndex(diag,l0));
        CHECK(hasIndex(diag,l2));
        }

    //SECTION("Bulk Case - IQTensor")
    //    {
    //    auto Op1 = randomITensor(QN(),S1,prime(S1),H0,H1);
    //    auto Op2 = randomITensor(QN(),S2,prime(S2),H1,H2);
    //    auto L = randomITensor(QN(),L0,prime(L0),H0);
    //    auto R = randomITensor(QN(),L2,prime(L2),H2);
    //    auto diag = lop.diag();
    //    CHECK(hasIndex(diag,S1));
    //    CHECK(hasIndex(diag,S2));
    //    CHECK(hasIndex(diag,L0));
    //    CHECK(hasIndex(diag,L2));
    //    }
    }
}


TEST_CASE("LocalMPO")
{
SECTION("LocalMPO As MPS")
    {
    auto N = 10;
    auto sites = SpinHalf(N,{"ConserveQNs=",true});

    auto ferro = InitState(sites,"Up");
    auto neel = InitState(sites);
    for(int j = 1; j <= N; ++j)
        {
        neel.set(j,j%2==1 ? "Up" : "Dn");
        }

    auto psiF = MPS(ferro);
    auto psiN = MPS(neel);

    auto lmps = LocalMPO(psiN);

    lmps.position(3,psiF);

    lmps.numCenter(1);
    lmps.position(6,psiF);
    
    lmps.numCenter(0);
    lmps.position(2,psiF);
    }
}


