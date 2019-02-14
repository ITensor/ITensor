#include "test.h"
#include "itensor/mps/mps.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinless.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
using std::vector;

TEST_CASE("MPSTest")
{

static const int N = 10;
SpinHalf shsites(N,{"ConserveQNs=",false});
SpinHalf shsitesQNs(N,{"ConserveQNs=",true});

InitState shFerro(shsites,"Up");
InitState shFerroQNs(shsitesQNs,"Up");
InitState shNeel(shsites);

for(int j = 1; j <= N; ++j)
    {
    shNeel.set(j,j%2==1 ? "Up" : "Dn");
    }

SECTION("Constructors")
    {
    }


//SECTION("QNCheck")
//    {
//    IQMPS psiNeel(shNeel);
//    CHECK(checkQNs(psiNeel));
//
//    CHECK_EQUAL(totalQN(psiNeel),QN(0));
//
//    IQMPS psiFerro(shFerro);
//    CHECK(checkQNs(psiFerro));
//
//    CHECK_EQUAL(totalQN(psiFerro),QN(10));
//    }

SECTION("hasQNs")
    {
    CHECK(hasQNs(shFerroQNs));
    CHECK(not hasQNs(shFerro));
    }

SECTION("Constructors (m==1)")
    {
    auto psi = MPS(shsitesQNs);
    auto l2 = commonIndex(psi.A(2),psi.A(3));
    CHECK(1==l2.dim());
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsitesQNs(n)==findIndex(psi.A(n),format("n=%d",n)));
        }
    }

SECTION("Constructors (m>1)")
    {
    auto m = 4;
    auto psi = MPS(shsites,m);
    auto l2 = commonIndex(psi.A(2),psi.A(3));
    CHECK(m==l2.dim());

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    for(int n = 1; n <= N; ++n)
      randomize(psi.Aref(n));

    psi.position(1);
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    psi.position(N);
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }
    }

SECTION("Random constructors (m==1)")
    {
    auto psi = randomMPS(shsites);
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    psi.position(1);
    auto normpsi = norm(psi);
    CHECK(normpsi>0);
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    psi.position(N);
    CHECK_CLOSE(normpsi,norm(psi));
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shsites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }
    }

SECTION("Random constructors, QN conserved (m==1)")
    {
    auto psi = randomMPS(shFerroQNs);
    CHECK(norm(psi)>0);
    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(shFerroQNs(n)==findIndex(psi.A(n),format("n=%d",n)));
        }
    }

SECTION("MPSAddition 1")
    {
    Spinless sites(10,{"ConserveQNs=",true});

    InitState i1(sites,"Emp"),
              i2(sites,"Emp");

    i1.set(1,"Occ");
    i2.set(2,"Occ");

    //"Valence bond" between sites 1 and 2
    MPS psi = ISqrt2*sum(MPS(i1),MPS(i2));

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    CHECK_CLOSE(norm(psi),1);
    CHECK_EQUAL(totalQN(psi),QN({"Nf",1}));
    }

SECTION("MPSAddition 2")
    {
    auto sites = SiteSet(10,2);
    auto psi1 = randomMPS(sites);
    auto psi2 = randomMPS(sites);
    auto psi = sum(psi1,psi2);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    CHECK_EQUAL(rank(psi.A(1)),2);
    CHECK_EQUAL(rank(psi.A(2)),3);
    CHECK_EQUAL(rank(psi.A(5)),3);
    CHECK_EQUAL(rank(psi.A(9)),3);
    CHECK_EQUAL(rank(psi.A(10)),2);
    }

//SECTION("PositionTest")
//    {
//    Spinless sites(10);
//
//    InitState init(sites,"Emp");
//    init.set(2,"Occ");
//    init.set(4,"Occ");
//    init.set(6,"Occ");
//
//    IQMPS psi(init);
//    psi.Anc(1) *= Complex_i;
//
//    psi.position(1,"Cutoff=1E-8");
//    CHECK_EQUAL(findCenter(psi),1);
//
//    psi.position(4,"Cutoff=1E-8");
//    CHECK_EQUAL(findCenter(psi),4);
//    }

SECTION("Orthogonalize")
    {
    auto N = 10;
    auto m = 20;
    auto sites = SpinHalf(10,{"ConserveQNs=",false});
    auto psi = MPS(sites);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    //Make a random MPS of bond dim. m
    auto links = vector<Index>(N+1);
    for(auto n : range1(N))
        {
        links.at(n) = Index(m,format("Link,l=%d",n));
        }
    psi.Aref(1) = randomITensor(links.at(1),sites(1));
    for(auto n : range1(2,N-1))
        {
        psi.Aref(n) = randomITensor(links.at(n-1),sites(n),links.at(n));
        }
    psi.Aref(N) = randomITensor(links.at(N-1),sites(N));

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    //Normalize psi
    auto n2 = overlap(psi,psi);
    psi.Aref(1) /= sqrt(n2);

    auto opsi = psi;

    //for(auto b : range1(psi.N()-1))
    //    {
    //    Print(linkInd(psi,b));
    //    }

    psi.orthogonalize({"Cutoff",1E-16});
    CHECK_CLOSE(overlap(opsi,psi),1.0);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    //for(auto b : range1(psi.N()-1))
    //    {
    //    Print(linkInd(psi,b));
    //    }

    for(int n = N; n > 1; --n)
        {
        auto li = commonIndex(psi.A(n),psi.A(n-1),"Link");
        auto rho = psi.A(n) * dag(prime(psi.A(n),li));
        auto id = ITensor(li,prime(li));
        for(auto l : range1(li.dim()))
            {
            id.set(li(l),prime(li)(l),1.0);
            }
        CHECK(norm(rho-id) < 1E-10);
        }

    psi.orthogonalize({"Maxm=",10,"Cutoff=",1E-16});
    for(auto b : range1(psi.N()-1))
        {
        CHECK(linkInd(psi,b).dim() <= 10);
        }

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    }

SECTION("Overlap - 1 site")
    {
    auto psi = MPS(1);
    auto s = Index(2,"s");
    psi.Aref(1) = randomITensor(s);
    CHECK_CLOSE(overlap(psi,psi),(psi.A(1)*psi.A(1)).real());
    }

}
