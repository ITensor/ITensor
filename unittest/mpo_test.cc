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
    auto sites = SpinHalf(10);
    auto W = MPO(sites);

    //Make a random MPS of bond dim. m
    auto links = vector<Index>(N+1);
    for(auto n : range1(N))
        {
        links.at(n) = Index(nameint("l",n),m);
        }
    W.Aref(1) = randomTensor(links.at(1),sites(1),prime(sites(1)));
    for(auto n : range1(2,N-1))
        {
        W.Aref(n) = randomTensor(links.at(n-1),sites(n),prime(sites(n)),links.at(n));
        }
    W.Aref(N) = randomTensor(links.at(N-1),sites(N),prime(sites(N)));

    //Normalize W
    auto n2 = overlap(W,W);
    W.Aref(1) /= sqrt(n2);

    auto oW = W;

    W.orthogonalize();

    CHECK_CLOSE(overlap(oW,W),1.0);

    for(int n = N; n > 1; --n)
        {
        auto li = commonIndex(W.A(n),W.A(n-1),Link);
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


    auto makeInds = [N](std::string name) -> vector<IQIndex>
        {
        auto ll = vector<IQIndex>(N);
        for(auto n : range1(N-1))
            {
            ll.at(n) = IQIndex(nameint(name,n),
                               Index("a",2),QN("Sz=",-1,"Nf=",-1),
                               Index("a",2),QN("Sz=",-1,"Nf=",+1),
                               Index("b",2),QN("Sz=",-1,"Nf=",0),
                               Index("c",2),QN("Sz=",+1,"Nf=",0),
                               Index("d",2),QN("Sz=",+1,"Nf=",-1),
                               Index("d",2),QN("Sz=",+1,"Nf=",+1));
            }
        return ll;
        };

    auto l1 = makeInds("I1_");
    auto l2 = makeInds("I2_");

    auto Z = QN("Sz=",0,"Nf=",0);

    auto A = IQMPO(sites);
    auto B = IQMPO(sites);
    A.Aref(1) = randomTensor(Z,sites(1),l1.at(1));
    B.Aref(1) = randomTensor(Z,sites(1),l2.at(1));
    for(int n = 2; n < N; ++n)
        {
        A.Aref(n) = randomTensor(Z,sites(n),dag(l1.at(n-1)),l1.at(n));
        B.Aref(n) = randomTensor(Z,sites(n),dag(l2.at(n-1)),l2.at(n));
        }
    A.Aref(N) = randomTensor(Z,sites(N),dag(l1.at(N-1)));
    B.Aref(N) = randomTensor(Z,sites(N),dag(l2.at(N-1)));

    auto C = sum(A,B);

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

    auto A = IQMPO(sites);
    auto Ia = IQIndex("I",Index("1",2),QN("Sz",1,"Nf",-1),
                          Index("1",1),QN("Sz",-1,"Nf",-1));
    A.Aref(1) = randomTensor(QN("Sz",-1,"Nf",1), prime(sites(1)), dag(Ia), dag(sites(1)));
    A.Aref(2) = randomTensor(QN("Sz",1,"Nf",-1), Ia, dag(sites(2)), prime(sites(2)));

    auto B = IQMPO(sites);
    auto Ib = IQIndex("I",Index("1",2),QN("Sz",1,"Nf",-1),
                          Index("1",1),QN("Sz",-1,"Nf",-1));
    B.Aref(1) = randomTensor(QN("Sz",0,"Nf",0), prime(sites(1)), dag(Ib), dag(sites(1)));
    B.Aref(2) = randomTensor(QN("Sz",0,"Nf",0), prime(sites(2)), Ib, dag(sites(2)));

    REQUIRE_NOTHROW(A.plusEq(B));
    }

SECTION("applyMPO (DensityMatrix)")
    {

    auto method = "DensityMatrix";

    auto N = 10;
    auto sites = SpinHalf(N);

    auto psi = MPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = MPO(ampo);
    auto K = MPO(ampo);
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

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000});

    CHECK_EQUAL(checkMPOProd(Hpsi,H,psi,1E-10),true);

    }

SECTION("applyMPO (Fit)")
    {

    auto method = "Fit";

    auto N = 10;
    auto sites = SpinHalf(N);

    auto psi = MPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = MPO(ampo);
    auto K = MPO(ampo);
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

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000,"Sweeps=",100});

    CHECK_EQUAL(checkMPOProd(Hpsi,H,psi,1E-10),true);

    // Now with a trial starting state
    auto Hpsi_2 = applyMPO(H,psi,Hpsi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000,"Sweeps=",100});

    CHECK_EQUAL(checkMPOProd(Hpsi_2,H,psi,1E-10),true);

    }

SECTION("errorMPOProd Scaling")
    {

    auto method = "DensityMatrix";

    auto N = 10;
    auto sites = SpinHalf(N);

    auto psi = MPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = MPO(ampo);
    auto K = MPO(ampo);
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

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"Maxm=",5000});

    //<Hpsi|Hpsi> is ~ 1E20, but normalization should take care of that
    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.);

    }

SECTION("Overlap <psi|HK|phi>")
    {
    detail::seed_quickran(1);

    auto N = 10;
    auto sites = SpinHalf(N);

    auto psi = MPS(sites);
    auto phi = MPS(sites);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = MPO(ampo);
    auto Hdag = H;
    auto K = H;
    //Randomize the MPOs to make sure they are non-Hermitian
    for(auto j : range1(N))
        {
        randomize(H.Aref(j));
        randomize(K.Aref(j));
        H.Aref(j) *= 0.2;
        K.Aref(j) *= 0.3;
        Hdag.Aref(j) = dag(swapPrime(H.A(j),0,1,Site));
        }

    auto Hdphi = applyMPO(Hdag,phi,{"Cutoff=",1E-13,"Maxm=",5000,"Method=","DensityMatrix"});
    auto Kpsi = applyMPO(K,psi,{"Cutoff=",1E-13,"Maxm=",5000,"Method=","DensityMatrix"});

    //Print(overlap(phi,H,K,psi));
    //Print(overlap(Hdphi,Kpsi));
    CHECK_CLOSE(overlap(phi,H,K,psi),overlap(Hdphi,Kpsi));
    }

SECTION("toMPO function")
    {
    auto N = 50;
    auto sites = Hubbard(N);

    auto makeInds = [N](std::string name) -> vector<IQIndex>
        {
        auto ll = vector<IQIndex>(N);
        for(auto n : range1(N-1))
            {
            ll.at(n) = IQIndex(nameint(name,n),
                               Index("a",2),QN("Sz=",-1,"Nf=",-1),
                               Index("a",2),QN("Sz=",-1,"Nf=",+1),
                               Index("b",2),QN("Sz=",-1,"Nf=",0),
                               Index("c",2),QN("Sz=",+1,"Nf=",0),
                               Index("d",2),QN("Sz=",+1,"Nf=",-1),
                               Index("d",2),QN("Sz=",+1,"Nf=",+1));
            }
        return ll;
        };

    auto ll = makeInds("I_");

    auto Z = QN("Sz=",0,"Nf=",0);

    auto A = IQMPO(sites);
    A.Aref(1) = randomTensor(Z,sites(1),ll.at(1));
    for(int n = 2; n < N; ++n)
        {
        A.Aref(n) = randomTensor(Z,sites(n),dag(ll.at(n-1)),ll.at(n));
        }
    A.Aref(N) = randomTensor(Z,sites(N),dag(ll.at(N-1)));

    auto a = toMPO(A);

    for(auto n : range1(N))
        {
        CHECK(norm(a.A(n) - ITensor(A.A(n))) < 1E-10);
        }

    }


SECTION("nmultMPO")
    {
    detail::seed_quickran(1);

    auto N = 4;
    auto sites = SpinHalf(N);

    //Use AutoMPO as a trick to get
    //an MPO with bond dimension > 1
    auto ampo = AutoMPO(sites);
    for(auto j : range1(N-1))
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    auto H = MPO(ampo);
    auto K = MPO(ampo);
    //Randomize the MPOs to make sure they are non-Hermitian
    for(auto j : range1(N))
        {
        randomize(H.Aref(j));
        randomize(K.Aref(j));
        H.Aref(j) *= 0.2;
        K.Aref(j) *= 0.3;
        }

    auto nsites = SpinHalf(N);
    for(auto j : range1(N))
        {
        K.Aref(j) *= prime(delta(Index(nsites(j)),sites(j)));
        }

    MPO R;
    nmultMPO(H,K,R,{"Cutoff=",1E-10});

    auto randomMPS = [](SiteSet const& sites, 
              int m)
        {
        auto psi = MPS(sites);
        auto N = sites.N();
        auto links = std::vector<Index>(N+1);
        for(int j = 0; j <= N; ++j)
            {
            links.at(j) = Index(format("l_%d",j),m);
            }
        for(int j = 1; j <= N; ++j)
            {
            psi.Aref(j) = randomTensor(links.at(j-1),sites(j),links.at(j));
            psi.Aref(j) /= norm(psi.A(j));
            }
        psi.Aref(1) *= randomTensor(links.at(0));
        psi.Aref(N) *= randomTensor(links.at(N));
        //printfln("<rpsi|rpsi> = %.12f",overlap(psi,psi));
        psi.position(1);
        psi.Aref(1) /= norm(psi.A(1));
        return psi;
        };

    //Print(R);

    int Ntest = 4;
    for(int n = 1; n <= Ntest; ++n)
        {
        auto psi = randomMPS(nsites,4);
        auto phi = randomMPS(sites,4);
        auto oR = overlap(psi,R,phi);
        auto oKH = overlap(psi,K,H,phi);
        CHECK_DIFF(oR,oKH,1E-5);
        }
    }
SECTION("Random MPOs")
      {
      auto sites = SpinHalf(6);
      /* create an MPO with a random tensor on every site with dimension k on every bond
       * optional argument "Orthogonalize":true in order to 
       * automatically orthogonalize after creation
       */
      int k = 2;
      auto H1 = MPO(sites,k);

      auto OrthoH1 = MPO(sites,k,{"Orthogonalize",true}); 

      /* create an MPO with a random tensor on every site 
       * with each bond dimension specified by a vector
       */
      std::vector<int> bondList = {2,3,5,3,2};
      auto H2 = MPO(sites,bondList);
      //bonds directly from initializer list
      auto H3 = MPO(sites,{3,3,3,3,3},{"Orthogonalize",true,"Scale",0.2});
      
      CHECK_EQUAL(maxM(H1),2);
      CHECK_EQUAL(maxM(H2),5);
      CHECK_EQUAL(maxM(H3),3);
      CHECK_CLOSE(averageM(H1),2);
      CHECK_CLOSE(averageM(H2),3);
      CHECK_CLOSE(averageM(H3),3);
      REQUIRE_NOTHROW(overlap(OrthoH1,OrthoH1));
      REQUIRE_NOTHROW(overlap(H2,H2));
      REQUIRE_NOTHROW(overlap(H3,H3));
      }
}
