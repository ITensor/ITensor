#include "test.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/util/print_macro.h"
#include "itensor/mps/sites/hubbard.h"
#include "itensor/mps/autompo.h"

using namespace itensor;
using namespace std;

TEST_CASE("MPOTest")
{

SECTION("Orthogonalize")
    {
    auto N = 10;
    auto m = 4;
    auto sites = SpinHalf(10,{"ConserveQNs=",false});
    auto W = MPO(sites);

    //Make a random MPO of bond dim. m
    auto links = vector<Index>(N+1);
    for(auto n : range1(N))
        {
        links.at(n) = Index(m,format("Link,l=%d",n));
        }
    W.Aref(1) = randomITensor(links.at(1),sites(1),prime(sites(1)));
    for(auto n : range1(2,N-1))
        {
        W.Aref(n) = randomITensor(links.at(n-1),sites(n),prime(sites(n)),links.at(n));
        }
    W.Aref(N) = randomITensor(links.at(N-1),sites(N),prime(sites(N)));

    //Normalize W
    auto n2 = overlap(W,W);
    W.Aref(1) /= sqrt(n2);

    auto oW = W;

    W.orthogonalize();

    CHECK_CLOSE(overlap(oW,W),1.0);

    for(int n = N; n > 1; --n)
        {
        auto li = commonIndex(W.A(n),W.A(n-1),"Link");
        CHECK(li==findIndex(W.A(n),format("l=%d",n-1)));
        CHECK(li==findIndex(W.A(n-1),format("l=%d",n-1)));
        CHECK(sites(n)==findIndex(W.A(n),format("n=%d",n),0));
        auto rho = W.A(n) * dag(prime(W.A(n),li));
        auto id = ITensor(li,prime(li));
        for(auto l : range1(li.m()))
            {
            id.set(li(l),prime(li)(l),1.0);
            }
        CHECK(norm(rho-id) < 1E-10);
        }
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
    A.Aref(1) = randomITensor(Z,sites(1),l1.at(1));
    B.Aref(1) = randomITensor(Z,sites(1),l2.at(1));
    for(int n = 2; n < N; ++n)
        {
        A.Aref(n) = randomITensor(Z,sites(n),dag(l1.at(n-1)),l1.at(n));
        B.Aref(n) = randomITensor(Z,sites(n),dag(l2.at(n-1)),l2.at(n));
        }
    A.Aref(N) = randomITensor(Z,sites(N),dag(l1.at(N-1)));
    B.Aref(N) = randomITensor(Z,sites(N),dag(l2.at(N-1)));

    auto C = sum(A,B);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(C.A(n),C.A(n+1),"Link");
        CHECK(ln==findIndex(C.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(C.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(C.A(n),format("n=%d",n),0));
        }

    auto AA = overlap(A,A);
    auto AB = overlap(A,B);
    auto AC = overlap(A,C);
    auto BB = overlap(B,B);
    auto BC = overlap(B,C);
    auto CC = overlap(C,C);

    // |(A+B)-C|^2 = (A+B-C)*(A+B-C) = A*A+2A*B-2A*C+B*B-2B*C+C*C

    auto diff2 = AA+2*AB-2*AC+BB-2*BC+CC;
    CHECK(diff2 < 1E-12);
    }

SECTION("Regression Test")
    {
    auto sites = Hubbard(2);

    auto A = MPO(sites);
    auto Ia = Index(QN({"Sz",1},{"Nf",-1,-1}),2,
                    QN({"Sz",-1},{"Nf",-1,-1}),1,"I");
    A.Aref(1) = randomITensor(QN({"Sz",-1},{"Nf",1,-1}), prime(sites(1)), dag(Ia), dag(sites(1)));
    A.Aref(2) = randomITensor(QN({"Sz",1},{"Nf",-1,}), Ia, dag(sites(2)), prime(sites(2)));

    auto B = MPO(sites);
    auto Ib = Index(QN({"Sz",1},{"Nf",-1,-1}),2,
                    QN({"Sz",-1},{"Nf",-1,-1}),1,"I");
    B.Aref(1) = randomITensor(QN({"Sz",0},{"Nf",0,-1}), prime(sites(1)), dag(Ib), dag(sites(1)));
    B.Aref(2) = randomITensor(QN({"Sz",0},{"Nf",0,-1}), prime(sites(2)), Ib, dag(sites(2)));

    REQUIRE_NOTHROW(A.plusEq(B));
    }

SECTION("applyMPO (DensityMatrix)")
    {

    auto method = "DensityMatrix";

    auto N = 10;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto psi = randomMPS(sites);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n)));
        }

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = toMPO(ampo);
    auto K = toMPO(ampo);

    //Randomize the MPOs to make sure they are non-Hermitian
    for(auto j : range1(N))
        {
        randomize(H.Aref(j));
        randomize(K.Aref(j));
        H.Aref(j) *= 0.2;
        K.Aref(j) *= 0.2;
        }

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(H.A(n),H.A(n+1),"Link");
        CHECK(ln==findIndex(H.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(H.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(H.A(n),format("n=%d",n),0));
        }

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"Maxm=",100});
    psi /= norm(psi);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n),0));
        }

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000});

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(Hpsi.A(n),Hpsi.A(n+1),"Link");
        CHECK(ln==findIndex(Hpsi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(Hpsi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(Hpsi.A(n),format("n=%d",n),0));
        }

    CHECK_EQUAL(checkMPOProd(Hpsi,H,psi,1E-10),true);

    }

SECTION("applyMPO (Fit)")
    {

    auto method = "Fit";

    auto N = 10;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto psi = randomMPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = toMPO(ampo);
    auto K = toMPO(ampo);
    //Randomize the MPOs to make sure they are non-Hermitian
    for(auto j : range1(N))
        {
        randomize(H.Aref(j));
        randomize(K.Aref(j));
        H.Aref(j) *= 0.2;
        K.Aref(j) *= 0.2;
        }

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"Maxm=",100});
    psi /= norm(psi);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n),0));
        }

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000,"Sweeps=",100});

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(Hpsi.A(n),Hpsi.A(n+1),"Link");
        CHECK(ln==findIndex(Hpsi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(Hpsi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(Hpsi.A(n),format("n=%d",n),0));
        }

    CHECK_EQUAL(checkMPOProd(Hpsi,H,psi,1E-10),true);

    // Now with a trial starting state
    auto Hpsi_2 = applyMPO(H,psi,Hpsi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000,"Sweeps=",100});

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(Hpsi_2.A(n),Hpsi_2.A(n+1),"Link");
        CHECK(ln==findIndex(Hpsi_2.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(Hpsi_2.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(Hpsi_2.A(n),format("n=%d",n),0));
        }

    CHECK_EQUAL(checkMPOProd(Hpsi_2,H,psi,1E-10),true);

    }

SECTION("errorMPOProd Scaling")
    {

    auto method = "DensityMatrix";

    auto N = 10;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto psi = randomMPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = toMPO(ampo);
    auto K = toMPO(ampo);
    //Randomize the MPOs to make sure they are non-Hermitian
    for(auto j : range1(N))
        {
        randomize(H.Aref(j));
        randomize(K.Aref(j));
        H.Aref(j) *= 10.0; //crazy large tensor
        K.Aref(j) *= 10.0;
        }

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"Maxm=",100});
    psi /= norm(psi);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(psi.A(n),psi.A(n+1),"Link");
        CHECK(ln==findIndex(psi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(psi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(psi.A(n),format("n=%d",n),0));
        }

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000});

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(Hpsi.A(n),Hpsi.A(n+1),"Link");
        CHECK(ln==findIndex(Hpsi.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(Hpsi.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(Hpsi.A(n),format("n=%d",n),0));
        }

    //<Hpsi|Hpsi> is ~ 1E20, but normalization should take care of that
    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.);

    }

//TODO: test this without using applyMPO()?
SECTION("Overlap <psi|HK|phi>")
    {
    detail::seed_quickran(1);

    auto N = 10;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});

    auto psi = randomMPS(sites);
    auto phi = randomMPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = toMPO(ampo);
    auto Hdag = H;
    auto K = H;
    //Randomize the MPOs to make sure they are non-Hermitian
    for(auto j : range1(N))
        {
        randomize(H.Aref(j));
        randomize(K.Aref(j));
        H.Aref(j) *= 0.2;
        K.Aref(j) *= 0.3;
        Hdag.Aref(j) = dag(swapPrime(H.A(j),0,1,"Site"));
        }

    auto Hdphi = applyMPO(Hdag,phi,{"Cutoff=",1E-13,"Maxm=",5000,"Method=","DensityMatrix"});
    auto Kpsi = applyMPO(K,psi,{"Cutoff=",1E-13,"Maxm=",5000,"Method=","DensityMatrix"});

    //Print(overlap(phi,H,K,psi));
    //Print(overlap(Hdphi,Kpsi));
    CHECK_CLOSE(overlap(phi,H,K,psi),overlap(Hdphi,Kpsi));
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
    A.Aref(1) = randomITensor(Z,sites(1),ll.at(1));
    for(int n = 2; n < N; ++n)
        {
        A.Aref(n) = randomITensor(Z,sites(n),dag(ll.at(n-1)),ll.at(n));
        }
    A.Aref(N) = randomITensor(Z,sites(N),dag(ll.at(N-1)));

    auto a = removeQNs(A);

    for(int n = 1; n < N; ++n)
        {
        auto ln = commonIndex(a.A(n),a.A(n+1),"Link");
        CHECK(ln==findIndex(a.A(n),format("l=%d",n)));
        CHECK(ln==findIndex(a.A(n+1),format("l=%d",n)));
        CHECK(sites(n)==findIndex(a.A(n),format("n=%d",n),0));
        }

    for(auto n : range1(N))
        {
        CHECK(norm(a.A(n) - removeQNs(A.A(n))) < 1E-10);
        }

    }

}
