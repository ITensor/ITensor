#include "test.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/sites/spinone.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/str.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/autompo.h"
#include "itensor/mps/dmrg.h"
#include "mps_mpo_test_helper.h"

using namespace itensor;
using namespace std;

TEST_CASE("MPOTest")
{

SECTION("Orthogonalize")
    {
    auto N = 20;
    auto d = 4;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
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
    auto n2 = trace(dag(W),W);
    W.ref(1) /= sqrt(n2);

    auto oW = W;

    oW.orthogonalize();

    CHECK_CLOSE(trace(dag(oW),W),1.0);
    CHECK(checkOrtho(oW));
    }

SECTION("Add MPOs")
    {
    auto N = 50;
    auto sites = Electron(N);

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

    auto AA = traceC(dag(A),A);
    auto AB = traceC(dag(A),B);
    auto AC = traceC(dag(A),C);
    auto BB = traceC(dag(B),B);
    auto BC = traceC(dag(B),C);
    auto CC = traceC(dag(C),C);

    // |(A+B)-C|^2 = (A+B-C)*(A+B-C) = A*A+2A*B-2A*C+B*B-2B*C+C*C

    auto diff2 = AA+2*AB-2*AC+BB-2*BC+CC;
    CHECK(std::abs(diff2) < 1E-12);
    }

SECTION("Regression Test")
    {
    auto sites = Electron(2);

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

    auto N = 20;
    auto sites = SpinHalf(N);
    auto initstate = InitState(sites,"Up");
    for( auto j : range1(N) ) if( j%2 == 1 )
      initstate.set(j,"Dn");

    auto psi = randomMPS(initstate,{"Complex=",true});

    CHECK(checkTags(psi));

    auto H = randomUnitaryMPO(sites);
    auto K = randomUnitaryMPO(sites);

    CHECK(checkTags(H));
    CHECK(checkTags(K));

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",200});

    CHECK(checkTags(psi,"Site,1","Link,0"));

    psi.noPrime("Site");

    CHECK(checkTags(psi));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",200});

    CHECK(checkTags(Hpsi,"Site,1","Link,0"));

    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.0);
    }

SECTION("applyMPO (DensityMatrix) with custom tags")
    {
    auto method = "DensityMatrix";

    auto N = 20;
    auto sites = SpinHalf(N);
    auto initstate = InitState(sites,"Up");
    for( auto j : range1(N) ) if( j%2 == 1 )
      initstate.set(j,"Dn");

    auto psi = randomMPS(initstate,{"Complex=",true});

    CHECK(checkTags(psi));

    psi.replaceTags("Site,S=1/2","MySite,bra");
    psi.replaceTags("Link","MyLink,psi");

    CHECK(checkTags(psi,"MySite,bra","MyLink,psi"));

    auto H = randomUnitaryMPO(sites);

    CHECK(checkTags(H));

    H.replaceTags("Site,S=1/2,0","MySite,bra,0");
    H.replaceTags("Site,S=1/2,1","MySite,ket,0");
    H.replaceTags("Link","MyLink,H");

    CHECK(checkTags(H,"MySite,bra","MySite,ket","MyLink,H"));

    auto K = randomUnitaryMPO(sites);

    CHECK(checkTags(K));

    K.replaceTags("Site,S=1/2,0","MySite,bra,0");
    K.replaceTags("Site,S=1/2,1","MySite,ket,0");
    K.replaceTags("Link","MyLink,K");

    CHECK(checkTags(K,"MySite,bra","MySite,ket","MyLink,K"));

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",100});
    psi /= norm(psi);

    CHECK(checkTags(psi,"MySite,ket","MyLink,psi"));

    psi.replaceTags("ket","bra");

    CHECK(checkTags(psi,"MySite,bra","MyLink,psi"));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",200});

    CHECK(checkTags(Hpsi,"MySite,ket","MyLink,psi"));

    Hpsi.replaceTags("psi","Hpsi");

    CHECK(checkTags(Hpsi,"MySite,ket","MyLink,Hpsi"));

    auto psik = dag(psi);
    psik.replaceTags("bra","ket").addTags("ket","MyLink");

    CHECK(checkTags(psik,"MySite,ket","MyLink,psi,ket"));

    auto O1 = psi(1)*H(1)*psik(1);
    auto O2 = Hpsi(1)*psik(1);
    for( auto n : range1(2,N) )
      {
      O1 = O1*psi(n)*H(n)*psik(n);
      O2 = O2*Hpsi(n)*psik(n);
      }

    CHECK_CLOSE(norm(O1-O2),0);

    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.0);
    }

SECTION("applyMPO (Fit)")
    {
    auto method = "Fit";

    auto N = 20;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    for( auto j : range1(N) ) if( j%2 == 1 )
      initstate.set(j,"Dn");

    auto psi = randomMPS(initstate,{"Complex=",true});

    CHECK(checkTags(psi));

    auto H = randomUnitaryMPO(sites);
    auto K = randomUnitaryMPO(sites);

    CHECK(checkTags(H));
    CHECK(checkTags(K));

    // Apply K to psi to entangle psi
    auto maxdim = 200;
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",maxdim});
    psi /= norm(psi);

    CHECK(checkTags(psi,"Site,1","Link,0"));

    psi.noPrime("Site");

    CHECK(checkTags(psi));

    CHECK( maxLinkDim(psi) <= maxdim );

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",maxdim,"Nsweep=",10});

    CHECK(checkTags(Hpsi,"Site,1","Link,0"));
    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.0);

    CHECK( maxLinkDim(Hpsi) <= maxdim );

    // Now with a trial starting state
    Hpsi = applyMPO(H,psi,Hpsi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",maxdim,"Nsweep=",5});

    CHECK( maxLinkDim(Hpsi) <= maxdim );

    CHECK(checkTags(Hpsi,"Site,1","Link,0"));
    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.0);
    }

SECTION("errorMPOProd Scaling")
    {
    auto method = "DensityMatrix";

    auto N = 20;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    for( auto j : range1(N) ) if( j%2 == 1 )
      initstate.set(j,"Dn");

    auto psi = randomMPS(initstate,{"Complex=",true});
    auto H = randomUnitaryMPO(sites);
    auto K = randomUnitaryMPO(sites);

    // Scale the MPOs to make the norms very large
    for( auto j : range1(N) )
        {
        H.ref(j) *= 20.0;
        K.ref(j) *= 20.0;
        }

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",100});
    psi.noPrime("Site");
    psi /= norm(psi);

    CHECK(checkTags(psi));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",200});

    CHECK(checkTags(Hpsi,"Site,1","Link,0"));

    //<Hpsi|Hpsi> is ~ 1E20, but normalization should take care of that
    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.);
    }

SECTION("applyMPO (Fit) with custom tags")
    {
    auto method = "Fit";

    auto N = 20;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    for( auto j : range1(N) ) if( j%2 == 1 )
      initstate.set(j,"Dn");

    auto psi = randomMPS(initstate,{"Complex=",true});

    CHECK(checkTags(psi));

    psi.replaceTags("Site,S=1/2","MySite,bra");
    psi.replaceTags("Link","MyLink,psi");

    CHECK(checkTags(psi,"MySite,bra","MyLink,psi"));

    auto H = randomUnitaryMPO(sites);

    CHECK(checkTags(H));

    H.replaceTags("Site,S=1/2,0","MySite,bra,0");
    H.replaceTags("Site,S=1/2,1","MySite,ket,0");
    H.replaceTags("Link","MyLink,H");

    CHECK(checkTags(H,"MySite,bra","MySite,ket","MyLink,H"));

    auto K = randomUnitaryMPO(sites);

    CHECK(checkTags(K));

    K.replaceTags("Site,S=1/2,0","MySite,bra,0");
    K.replaceTags("Site,S=1/2,1","MySite,ket,0");
    K.replaceTags("Link","MyLink,K");

    CHECK(checkTags(K,"MySite,bra","MySite,ket","MyLink,K"));

    // Apply K to psi to entangle psi
    psi = applyMPO(K,psi,{"Cutoff=",0.,"MaxDim=",100});
    psi /= norm(psi);

    CHECK(checkTags(psi,"MySite,ket","MyLink,psi"));

    psi.replaceTags("ket","bra");

    CHECK(checkTags(psi,"MySite,bra","MyLink,psi"));

    auto Hpsi = applyMPO(H,psi,{"Method=",method,"Cutoff=",1E-13,"MaxDim=",200});

    CHECK(checkTags(Hpsi,"MySite,ket","MyLink,psi"));

    Hpsi.replaceTags("psi","Hpsi");

    CHECK(checkTags(Hpsi,"MySite,ket","MyLink,Hpsi"));

    auto psik = dag(psi);
    psik.replaceTags("bra","ket").addTags("ket","MyLink");

    CHECK(checkTags(psik,"MySite,ket","MyLink,psi,ket"));

    auto O1 = psi(1)*H(1)*psik(1);
    auto O2 = Hpsi(1)*psik(1);
    for( auto n : range1(2,N) )
      {
      O1 = O1*psi(n)*H(n)*psik(n);
      O2 = O2*Hpsi(n)*psik(n);
      }

    CHECK_CLOSE(norm(O1-O2),0);

    CHECK_CLOSE(errorMPOProd(Hpsi,H,psi),0.0);
    }

SECTION("Inner <Hpsi|Kphi> and <psi|H^{d}K|phi>")
    {
    detail::seed_quickran(1);

    auto N = 20;
    auto sites = SpinHalf(N);

    auto initstate = InitState(sites,"Up");
    for( auto j : range1(N) ) if( j%2 == 1 )
      initstate.set(j,"Dn");

    auto psi = randomMPS(initstate,{"Complex=",true});
    auto phi = randomMPS(initstate,{"Complex=",true});
    auto H = randomUnitaryMPO(sites);
    auto K = randomUnitaryMPO(sites);

    CHECK(checkTags(H));
    CHECK(checkTags(K));

    auto Hphi = applyMPO(H,phi,{"Cutoff=",1E-13,"MaxDim=",200,"Method=","DensityMatrix"});
    auto Kpsi = applyMPO(K,psi,{"Cutoff=",1E-13,"MaxDim=",200,"Method=","DensityMatrix"});

    CHECK(checkTags(Hphi,"Site,1","Link,0"));
    CHECK(checkTags(Kpsi,"Site,1","Link,0"));

    // Check inner(A,x,B,y) = <Ax|By>
    CHECK_CLOSE(innerC(H,phi,K,psi),innerC(Hphi,Kpsi));
    CHECK_CLOSE(innerC(Hphi,K,psi),innerC(Hphi,Kpsi));

    // Check inner(x,A,B,y) = <x|AB|y>
    CHECK_CLOSE(innerC(phi,dag(H),K,psi),innerC(Hphi,Kpsi));
    CHECK_CLOSE(innerC(phi,dag(H),Kpsi),innerC(Hphi,Kpsi));
    }

SECTION("Remove QNs from MPO")
    {
    auto N = 50;
    auto sites = Electron(N);

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

  auto A = randomUnitaryMPO(sites);
  auto B = randomUnitaryMPO(sites);

  CHECK(checkTags(A));
  CHECK(checkTags(B));

  // Check nmultMPO with all shared site indices
  // throws an error
  CHECK_THROWS_AS(nmultMPO(A,B),ITError);

  // By default, C-links get the tags of A
  auto C = nmultMPO(A,prime(B));

  CHECK(checkTags(C,"Site,0","Site,2","Link,0"));

  // Check the product by calculating the trace
  CHECK_CLOSE(traceC(A,prime(B)),traceC(C));

  // Check the product by calculating expectation values
  auto initstate = InitState(sites,"Up");
  for( auto j : range1(N) ) if( j%2 == 1 )
    initstate.set(j,"Dn");
  auto Ncheck = 20;
  for( [[maybe_unused]] auto n : range1(Ncheck) )
    {
    auto V = randomMPS(initstate,{"Complex=",true});
    CHECK_CLOSE(innerC(V,C,V),innerC(V,prime(B),A,V));
    }
  }

SECTION("nmultMPO (custom tags)")
  {
  auto N = 4;
  auto sites = SpinHalf(N);
  auto A = randomUnitaryMPO(sites);
  auto B = randomUnitaryMPO(sites);

  // Set up some custom tags
  A.replaceTags("Site,S=1/2,0","x,0");
  A.replaceTags("Site,S=1/2,1","y,0");
  A.replaceTags("Link","Alink");
  B.replaceTags("Site,S=1/2,0","y,0");
  B.replaceTags("Site,S=1/2,1","z,0");
  B.replaceTags("Link","Blink");

  CHECK(checkTags(A,"x","y","Alink"));
  CHECK(checkTags(B,"y","z","Blink"));

  auto C = nmultMPO(A,B);

  CHECK(checkTags(C,"x","z","Alink"));

  CHECK_CLOSE(traceC(A,B),traceC(C));
  }

SECTION("DMRG")
  {
  int N = 32;
  auto sites = SpinHalf(N,{"ConserveQNs=",false});
  auto psi0 = randomMPS(InitState(sites,"Up"));

  auto h = 0.5; // Critical point

  auto ampo = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
      {
      ampo += -1.0,"Sx",j,"Sx",j+1;
      ampo += -h,"Sz",j;
      }
  ampo += -h,"Sz",N;    
  auto H = toMPO(ampo);

  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 10,20,30;
  sweeps.cutoff() = 1E-12;
  auto [Energy,psi] = dmrg(H,psi0,sweeps,{"Silent",true});
  auto energy = Energy/N;
  (void)psi;

  // Exact energy for transverse field Ising model
  // with open boundary conditions
  auto Energy_exact = 1.0 - 1.0/sin(Pi/(2*(2*N+1)));
  auto energy_exact = Energy_exact/(4*N);
  CHECK_CLOSE((energy-energy_exact)/energy_exact,0.);
  }

SECTION("DMRG with QNs")
  {
  int N = 32;
  auto sites = SpinHalf(N,{"ConserveSz=",false,
                           "ConserveParity=",true});
  auto psi0 = randomMPS(InitState(sites,"Up"));

  auto h = 0.5; // Critical point

  auto ampo = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
      {
      ampo += -1.0,"Sx",j,"Sx",j+1;
      ampo += -h,"Sz",j;
      }
  ampo += -h,"Sz",N;
  auto H = toMPO(ampo);

  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 10,20,30;
  sweeps.cutoff() = 1E-12;
  auto [Energy,psi] = dmrg(H,psi0,sweeps,{"Silent",true});
  auto energy = Energy/N;
  (void)psi;

  // Exact energy for transverse field Ising model
  // with open boundary conditions
  auto Energy_exact = 1.0 - 1.0/sin(Pi/(2*(2*N+1)));
  auto energy_exact = Energy_exact/(4*N);
  CHECK_CLOSE((energy-energy_exact)/energy_exact,0.);
  }

}
