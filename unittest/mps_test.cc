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
SpinHalf shsites(N);

InitState shFerro(shsites,"Up");
InitState shNeel(shsites);

for(int j = 1; j <= N; ++j)
    {
    shNeel.set(j,j%2==1 ? "Up" : "Dn");
    }

SECTION("Constructors")
    {
    }

SECTION("QNCheck")
    {
    IQMPS psiNeel(shNeel);
    CHECK(checkQNs(psiNeel));

    CHECK_EQUAL(totalQN(psiNeel),QN(0));

    IQMPS psiFerro(shFerro);
    CHECK(checkQNs(psiFerro));

    CHECK_EQUAL(totalQN(psiFerro),QN(10));
    }

//SECTION("MPSAddition")
//    {
//    Spinless sites(10);
//
//    InitState i1(sites,"Emp"),
//              i2(sites,"Emp");
//
//    i1.set(1,"Occ");
//    i2.set(2,"Occ");
//
//    //"Valence bond" between sites 1 and 2
//    MPS psi = ISqrt2*sum(MPS(i1),MPS(i2));
//
//    CHECK_CLOSE(norm(psi),1);
//
//    IQMPS iqpsi = ISqrt2*sum(IQMPS(i1),IQMPS(i2));
//
//    CHECK_EQUAL(totalQN(iqpsi),QN(0,1));
//    }

SECTION("PositionTest")
    {
    Spinless sites(10);

    InitState init(sites,"Emp");
    init.set(2,"Occ");
    init.set(4,"Occ");
    init.set(6,"Occ");

    IQMPS psi(init);
    psi.Anc(1) *= Complex_i;

    psi.position(1,"Cutoff=1E-8");
    CHECK_EQUAL(findCenter(psi),1);

    psi.position(4,"Cutoff=1E-8");
    CHECK_EQUAL(findCenter(psi),4);
    }

SECTION("Orthogonalize")
    {
    auto N = 10;
    auto m = 20;
    auto sites = SpinHalf(10);
    auto psi = MPS(sites);

    //Make a random MPS of bond dim. m
    auto links = vector<Index>(N+1);
    for(auto n : range1(N))
        {
        links.at(n) = Index(nameint("l",n),m);
        }
    psi.Aref(1) = randomTensor(links.at(1),sites(1));
    for(auto n : range1(2,N-1))
        {
        psi.Aref(n) = randomTensor(links.at(n-1),sites(n),links.at(n));
        }
    psi.Aref(N) = randomTensor(links.at(N-1),sites(N));

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

    //for(auto b : range1(psi.N()-1))
    //    {
    //    Print(linkInd(psi,b));
    //    }

    for(int n = N; n > 1; --n)
        {
        auto li = commonIndex(psi.A(n),psi.A(n-1),Link);
        auto rho = psi.A(n) * dag(prime(psi.A(n),li));
        auto id = ITensor(li,prime(li));
        for(auto l : range1(li.m()))
            {
            id.set(li(l),prime(li)(l),1.0);
            }
        CHECK(norm(rho-id) < 1E-10);
        }

    psi.orthogonalize({"Maxm=",10,"Cutoff=",1E-16});
    for(auto b : range1(psi.N()-1))
        {
        CHECK(linkInd(psi,b).m() <= 10);
        }

    }

SECTION("Overlap - 1 site")
    {
    auto psi = MPS(1);
    auto s = Index("s",2);
    psi.Aref(1) = randomTensor(s);
    CHECK_CLOSE(overlap(psi,psi),(psi.A(1)*psi.A(1)).real());
    }

}
