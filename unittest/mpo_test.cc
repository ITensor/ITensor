#include "test.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"
#include "itensor/mps/sites/hubbard.h"
#include "itensor/mps/autompo.h"
#include "mps_mpo_test_helper.h"

using namespace itensor;
using namespace std;

TEST_CASE("MPOTest")
{

SECTION("Orthogonalize")
    {
    auto N = 10;
    auto d = 4;
    auto sites = SpinHalf(10,{"ConserveQNs=",false});
    auto W = MPO(N);

    //Make a random MPO of bond dim. m
    auto links = vector<Index>(N+1);
    for(auto n : range1(N))
        links[n] = Index(d,format("MyLink,l=%d",n));
    W.ref(1) = randomITensorC(links[1],sites(1),prime(sites(1)));
    for(auto n : range1(2,N-1))
        W.ref(n) = randomITensorC(links[n-1],sites(n),prime(sites(n)),links[n]);
    W.ref(N) = randomITensorC(links[N-1],sites(N),prime(sites(N)));
    W.replaceTags("Site,0","MySite,bra,0");
    W.replaceTags("Site,1","MySite,ket,0");

    CHECK(checkTags(W,"MySite,bra","MySite,ket","MyLink"));

    //Normalize W
    auto n2 = overlap(W,W);
    W.ref(1) /= sqrt(n2);

    auto oW = W;

    oW.orthogonalize();

    CHECK_CLOSE(overlap(oW,W),1.0);
    CHECK(checkOrtho(oW));
    }

SECTION("Add MPOs")
    {
    auto N = 50;
    auto sites = Hubbard(N);

    auto makeInds = [N](std::string name) -> vector<Index>
        {
        auto ll = vector<Index>(N);
        for(auto n : range1(N-1))
            {
            auto ts = format("%s,l=%d",name,n);
            ll.at(n) = Index(QN({"Sz",-1},{"Nf",-1,-1}),2,
                             QN({"Sz",-1},{"Nf",+1,-1}),2,
                             QN({"Sz",-1},{"Nf=",0,-1}),2,
                             QN({"Sz",+1},{"Nf=",0,-1}),2,
                             QN({"Sz",+1},{"Nf=",-1,-1}),2,
                             QN({"Sz",+1},{"Nf=",+1,-1}),2,
                             ts);
            }
        return ll;
        };

    auto l1 = makeInds("Link,I1");
    auto l2 = makeInds("Link,I2");

    auto Z = QN({"Sz",0},{"Nf",0,-1});

    auto A = MPO(sites);
    auto B = MPO(sites);
    A.ref(1) = randomITensorC(Z,sites(1),prime(dag(sites(1))),l1.at(1));
    A.ref(1) /= norm(A(1));
    B.ref(1) = randomITensorC(Z,sites(1),prime(dag(sites(1))),l2.at(1));
    B.ref(1) /= norm(B(1));
    for(int n = 2; n < N; ++n)
        {
        A.ref(n) = randomITensorC(Z,sites(n),prime(dag(sites(n))),dag(l1.at(n-1)),l1.at(n));
        A.ref(n) /= norm(A(n));
        B.ref(n) = randomITensorC(Z,sites(n),prime(dag(sites(n))),dag(l2.at(n-1)),l2.at(n));
        B.ref(n) /= norm(B(n));
        }
    A.ref(N) = randomITensorC(Z,sites(N),prime(dag(sites(N))),dag(l1.at(N-1)));
    A.ref(N) /= norm(A(N));
    B.ref(N) = randomITensorC(Z,sites(N),prime(dag(sites(N))),dag(l2.at(N-1)));
    B.ref(N) /= norm(B(N));

    auto C = sum(A,B);

    // Check C gets the tags of A
    CHECK(checkTags(C,"Site,0","Site,1","Link,I1"));

    auto AA = overlapC(A,A);
    auto AB = overlapC(A,B);
    auto AC = overlapC(A,C);
    auto BB = overlapC(B,B);
    auto BC = overlapC(B,C);
    auto CC = overlapC(C,C);

    // |(A+B)-C|^2 = (A+B-C)*(A+B-C) = A*A+2A*B-2A*C+B*B-2B*C+C*C

    auto diff2 = AA+2*AB-2*AC+BB-2*BC+CC;
    CHECK(std::abs(diff2) < 1E-12);
    }

SECTION("Regression Test")
    {
    auto sites = Hubbard(2);

    auto A = MPO(sites);
    auto Ia = Index(QN({"Sz",1},{"Nf",-1,-1}),2,
                    QN({"Sz",-1},{"Nf",-1,-1}),1,"I");
    A.ref(1) = randomITensor(QN({"Sz",-1},{"Nf",1,-1}), prime(sites(1)), dag(Ia), dag(sites(1)));
    A.ref(2) = randomITensor(QN({"Sz",1},{"Nf",-1,}), Ia, dag(sites(2)), prime(sites(2)));

    auto B = MPO(sites);
    auto Ib = Index(QN({"Sz",1},{"Nf",-1,-1}),2,
                    QN({"Sz",-1},{"Nf",-1,-1}),1,"I");
    B.ref(1) = randomITensor(QN({"Sz",0},{"Nf",0,-1}), prime(sites(1)), dag(Ib), dag(sites(1)));
    B.ref(2) = randomITensor(QN({"Sz",0},{"Nf",0,-1}), prime(sites(2)), Ib, dag(sites(2)));

    REQUIRE_NOTHROW(A.plusEq(B));
    }

SECTION("applyMPO (DensityMatrix)")
    {
    auto method = "DensityMatrix";

    auto N = 10;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    auto psi = randomMPS(initstate,{"Complex=",true});

    CHECK(checkTags(psi));

    auto H = randomMPO(sites,{"Complex=",true});
    auto K = randomMPO(sites,{"Complex=",true});

    CHECK(checkTags(H));
    CHECK(checkTags(K));

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",100});
    psi /= norm(psi);

    CHECK(checkTags(psi));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",5000});

    CHECK(checkTags(Hpsi));
    CHECK_EQUAL(checkMPOProd(Hpsi,H,psi,1E-10),true);
    }

SECTION("applyMPO (Fit)")
    {
    auto method = "Fit";

    auto N = 10;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    auto psi = randomMPS(initstate,{"Complex=",true});

    CHECK(checkTags(psi));

    auto H = randomMPO(sites,{"Complex=",true});
    auto K = randomMPO(sites,{"Complex=",true});

    CHECK(checkTags(H));
    CHECK(checkTags(K));

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",100});
    psi /= norm(psi);

    CHECK(checkTags(psi));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",5000,"Nsweep=",100});

    CHECK(checkTags(Hpsi));
    CHECK(checkMPOProd(Hpsi,H,psi,1E-10));

    // Now with a trial starting state
    Hpsi = applyMPO(H,psi,Hpsi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",5000,"Nsweep=",100});

    CHECK(checkTags(Hpsi));
    CHECK(checkMPOProd(Hpsi,H,psi,1E-10));
    }

SECTION("errorMPOProd Scaling")
    {
    auto method = "DensityMatrix";

    auto N = 10;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    auto psi = randomMPS(initstate,{"Complex=",true});
    auto H = randomMPO(sites,{"Complex=",true});
    auto K = randomMPO(sites,{"Complex=",true});

    // Scale the MPOs to make the norms very large
    for( auto j : range1(N) )
        {
        H.ref(j) *= 20.0;
        K.ref(j) *= 20.0;
        }

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",100});
    psi /= norm(psi);

    CHECK(checkTags(psi));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",5000});

    CHECK(checkTags(Hpsi));

    //<Hpsi|Hpsi> is ~ 1E20, but normalization should take care of that
    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.);
    }

SECTION("Overlap <psi|HK|phi>")
    {
    detail::seed_quickran(1);

    auto N = 10;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    auto psi = randomMPS(initstate,{"Complex=",true});
    auto phi = randomMPS(initstate,{"Complex=",true});
    auto H = randomMPO(sites,{"Complex=",true});
    auto K = randomMPO(sites,{"Complex=",true});

    CHECK(checkTags(H));
    CHECK(checkTags(K));

    auto Hd = H;
    for(auto j : range1(N))
        {
        auto s = siteInds(Hd,j);
        Hd.ref(j) = dag(swapInds(H(j),{s(1)},{s(2)}));
        }

    auto Hdphi = applyMPO(Hd,phi,{"Cutoff=",1E-13,"MaxDim=",5000,"Method=","DensityMatrix"});
    auto Kpsi = applyMPO(K,psi,{"Cutoff=",1E-13,"MaxDim=",5000,"Method=","DensityMatrix"});

    CHECK(checkTags(Hdphi));
    CHECK(checkTags(Kpsi));
    CHECK_CLOSE(overlapC(phi,H,K,psi),overlapC(Hdphi,Kpsi));
    }

SECTION("Remove QNs from MPO")
    {
    auto N = 50;
    auto sites = Hubbard(N);

    auto makeInds = [N](std::string name) -> vector<Index>
        {
        auto ll = vector<Index>(N);
        for(auto n : range1(N-1))
            {
            auto ts = format("%s,l=%d",name,n);
            ll.at(n) = Index(QN({"Sz",-1},{"Nf",-1,-1}),2,
                             QN({"Sz",-1},{"Nf",+1,-1}),2,
                             QN({"Sz",-1},{"Nf",0,-1}),2,
                             QN({"Sz",+1},{"Nf",0,-1}),2,
                             QN({"Sz",+1},{"Nf",-1,-1}),2,
                             QN({"Sz",+1},{"Nf",+1,-1}),2,
                             ts);
            }
        return ll;
        };

    auto ll = makeInds("Link,I");

    auto Z = QN({"Sz",0},{"Nf",0,-1});

    auto A = MPO(sites);
    A.ref(1) = randomITensorC(Z,sites(1),prime(dag(sites(1))),ll.at(1));
    for(int n = 2; n < N; ++n)
        A.ref(n) = randomITensorC(Z,sites(n),prime(dag(sites(n))),dag(ll.at(n-1)),ll.at(n));
    A.ref(N) = randomITensorC(Z,sites(N),prime(dag(sites(N))),dag(ll.at(N-1)));

    CHECK(checkTags(A,"Site,0","Site,1","Link,I"));

    auto a = removeQNs(A);

    CHECK(checkTags(a,"Site,0","Site,1","Link,I"));
    for(auto n : range1(N))
        CHECK(norm(a(n) - removeQNs(A(n))) < 1E-10);
    }

SECTION("nmultMPO")
  {
  auto N = 4;
  auto sites = SpinHalf(N);

  auto A = randomMPO(sites,{"Complex=",true});
  auto B = randomMPO(sites,{"Complex=",true});

  // By default, C-links get the tags of A
  auto C = nmultMPO(A,B);

  // Check the product by calculating expectation values
  auto initstate = InitState(sites,"Up");
  auto V = randomMPS(initstate,{"Complex=",true});

  CHECK_CLOSE(overlapC(V,C,V),overlapC(V,B,A,V));
  }

//SECTION("nmultMPO (custom tags)")
//  {
//  auto N = 4;
//  auto sites = SpinHalf(N);
//  auto A = randomMPO(sites,{"Complex=",true});
//  auto B = randomMPO(sites,{"Complex=",true});
//
//  // Set up some custom tags
//  A.replaceTags("Site,S=1/2,0","x,0");
//  A.replaceTags("Site,S=1/2,1","y,0");
//  A.replaceTags("Link","Alink");
//  B.replaceTags("Site,S=1/2,0","y,0");
//  B.replaceTags("Site,S=1/2,1","z,0");
//  B.replaceTags("Link","Blink");
//
//  CHECK(checkTags(A,"x","y","Alink"));
//  CHECK(checkTags(B,"y","z","Blink"));
//
//  auto C = nmultMPO(A,B);
//
//  CHECK(checkTags(C,"x","z","Alink"));
//  CHECK(checkTags(Cr,"x","z","Blink"));
//
//  auto initstate = InitState(sites,"Up");
//  auto V = randomMPS(initstate);
//  V.replaceTags("Link","Vlink");
//  V.replaceTags("Site,S=1/2","x");
//
//  CHECK(checkTags(V,"x","Vlink"));
//
//  auto Vc = V;
//  Vc.dag().replaceTags("x","z").addTags("dag","Vlink");
//
//  CHECK(checkTags(Vc,"z","Vlink,dag"));
//
//  auto VcCV = V(1)*C(1)*Vc(1);
//  for( auto i : range1(2,N) ) { VcCV *= V(i)*C(i)*Vc(i); }
//  auto VcABV = V(1)*A(1)*B(1)*Vc(1);
//  for( auto i : range1(2,N) ) { VcABV *= V(i)*A(i)*B(i)*Vc(i); }
//
//  CHECK_CLOSE(eltC(VcCV),eltC(VcABV));
//  }

}
