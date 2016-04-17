#include "test.h"
#include "itensor/mps/localop.h"
#include "itensor/mps/localmpo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

TEST_CASE("LocalOp")
{
auto s1 = Index("s1",2,Site);
auto s2 = Index("s2",2,Site);
auto h0 = Index("h0",4);
auto h1 = Index("h1",4);
auto h2 = Index("h2",4);
auto l0 = Index("l0",10);
auto l2 = Index("l2",10);

auto S1 = IQIndex("S1",Index("s1-",1,Site),QN(-1),
                       Index("s1+",1,Site),QN(+1));
auto S2 = IQIndex("S2",Index("s2-",1,Site),QN(-1),
                       Index("s2+",1,Site),QN(+1));
auto H0 = IQIndex("H0",Index("h0-2",4),QN(-2),
                       Index("h0+0",8),QN(+0),
                       Index("h0+2",4),QN(+2));
auto H1 = IQIndex("H1",Index("h1-2",4),QN(-2),
                       Index("h1+0",8),QN(+0),
                       Index("h1+2",4),QN(+2));
auto H2 = IQIndex("H2",Index("h2-2",4),QN(-2),
                       Index("h2+0",8),QN(+0),
                       Index("h2+2",4),QN(+2));
auto L0 = IQIndex("L0",Index("l0-2",4),QN(-2),
                       Index("l0+0",8),QN(+0),
                       Index("l0+2",4),QN(+2));
auto L2 = IQIndex("L2",Index("l2-2",4),QN(-2),
                       Index("l2+0",8),QN(+0),
                       Index("l2+2",4),QN(+2));


SECTION("Product")
    {
    SECTION("Bulk Case")
        {
        auto Op1 = randomTensor(s1,prime(s1),h0,h1);
        auto Op2 = randomTensor(s2,prime(s2),h1,h2);
        auto L = randomTensor(l0,prime(l0),h0);
        auto R = randomTensor(l2,prime(l2),h2);
        auto lop = LocalOp<ITensor>(Op1,Op2,L,R);
        auto psi = randomTensor(l0,s1,s2,l2);
        auto Hpsi = ITensor();
        lop.product(psi,Hpsi);
        CHECK(hasindex(Hpsi,s1));
        CHECK(hasindex(Hpsi,s2));
        CHECK(hasindex(Hpsi,l0));
        CHECK(hasindex(Hpsi,l2));
        }
    }

SECTION("Diag")
    {
    SECTION("Bulk Case - ITensor")
        {
        auto Op1 = randomTensor(s1,prime(s1),h0,h1);
        auto Op2 = randomTensor(s2,prime(s2),h1,h2);
        auto L = randomTensor(l0,prime(l0),h0);
        auto R = randomTensor(l2,prime(l2),h2);
        auto lop = LocalOp<ITensor>(Op1,Op2,L,R);
        auto diag = lop.diag();
        CHECK(hasindex(diag,s1));
        CHECK(hasindex(diag,s2));
        CHECK(hasindex(diag,l0));
        CHECK(hasindex(diag,l2));
        }

    //SECTION("Bulk Case - IQTensor")
    //    {
    //    auto Op1 = randomTensor(QN(),S1,prime(S1),H0,H1);
    //    auto Op2 = randomTensor(QN(),S2,prime(S2),H1,H2);
    //    auto L = randomTensor(QN(),L0,prime(L0),H0);
    //    auto R = randomTensor(QN(),L2,prime(L2),H2);
    //    auto lop = LocalOp<IQTensor>(Op1,Op2,L,R);
    //    auto diag = lop.diag();
    //    CHECK(hasindex(diag,S1));
    //    CHECK(hasindex(diag,S2));
    //    CHECK(hasindex(diag,L0));
    //    CHECK(hasindex(diag,L2));
    //    }
    }
}


TEST_CASE("LocalMPO")
{
SECTION("LocalMPO As MPS")
    {
    auto N = 10;
    auto sites = SpinHalf(N);

    auto ferro = InitState(sites,"Up");
    auto neel = InitState(sites);
    for(int j = 1; j <= N; ++j)
        {
        neel.set(j,j%2==1 ? "Up" : "Dn");
        }

    auto psiF = IQMPS(ferro);
    auto psiN = IQMPS(neel);

    auto lmps = LocalMPO<IQTensor>(psiN);
    lmps.position(3,psiF);
    }
}


