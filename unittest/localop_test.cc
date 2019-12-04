#include "test.h"
#include "itensor/mps/localop.h"
#include "itensor/mps/localmpo.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/autompo.h"
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
        auto lop = LocalOp(Op1,Op2,L,R,{"NumCenter=",2});
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
        auto lop = LocalOp(Op1,L,R,{"NumCenter=",1});
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
        auto lop = LocalOp(L,R,{"NumCenter=",0});
        auto psi = randomITensor(l0,l1);
        auto Hpsi = ITensor();
        lop.product(psi,Hpsi);
        CHECK(hasIndex(Hpsi,l0));
        CHECK(hasIndex(Hpsi,l1));
        }

    }

SECTION("LocalOp: check product")
  {
  auto s1 = Index(2,"Site");
  auto s2 = Index(2,"Site");
  auto h0 = Index(2,"Link");
  auto h1 = Index(2,"Link");
  auto h2 = Index(2,"Link");
  auto l0 = Index(2,"Link");
  auto l1 = Index(2,"Link");
  auto l2 = Index(2,"Link");
  auto r0 = Index(2,"Link");
  auto r1 = Index(2,"Link");
  auto r2 = Index(2,"Link");


  auto Op1 = randomITensor(s1,prime(s1),h0,h1);
  auto Op2 = randomITensor(s2,prime(s2),h1,h2);
  auto L0 = randomITensor(l0,prime(l0),h0);
  auto L1 = randomITensor(l1,prime(l1),h1);
  auto L2 = randomITensor(l2,prime(l2),h2);
  auto R0 = randomITensor(r0,prime(r0),h0);
  auto R1 = randomITensor(r1,prime(r1),h1);
  auto R2 = randomITensor(r2,prime(r2),h2);

  auto psi0 = randomITensor(l0,r0);
  auto psi1 = randomITensor(l0,s1,r1);
  auto psi2 = randomITensor(l0,s1,s2,r2);

  auto H0 = LocalOp(L0,R0,{"NumCenter=",0});
  auto H1 = LocalOp(Op1,L0,R1,{"NumCenter=",1});
  auto H2 = LocalOp(Op1,Op2,L0,R2,{"NumCenter=",2});

  auto Hpsi0 = psi0;
  H0.product(psi0,Hpsi0);

  auto Hpsi1 = psi1;
  H1.product(psi1,Hpsi1);

  auto Hpsi2 = psi2;
  H2.product(psi2,Hpsi2);

  CHECK_CLOSE(norm(Hpsi0-noPrime(psi0*L0*R0)),0.);
  CHECK_CLOSE(norm(Hpsi1-noPrime(psi1*L0*Op1*R1)),0.);
  CHECK_CLOSE(norm(Hpsi2-noPrime(psi2*L0*Op1*Op2*R2)),0.);
  }

SECTION("Diag")
    {
    SECTION("Bulk Case - ITensor")
        {
        auto Op1 = randomITensor(s1,prime(s1),h0,h1);
        auto Op2 = randomITensor(s2,prime(s2),h1,h2);
        auto L = randomITensor(l0,prime(l0),h0);
        auto R = randomITensor(l2,prime(l2),h2);
        auto lop = LocalOp(Op1,Op2,L,R,{"NumCenter=",2});
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

SECTION("Test product")
  {
  int N = 10;
  auto sites = SpinHalf(N);
  auto ampo = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
      {
      ampo += 0.5,"S+",j,"S-",j+1;
      ampo += 0.5,"S-",j,"S+",j+1;
      ampo +=     "Sz",j,"Sz",j+1;
      }
  auto H = toMPO(ampo);
  auto psi = MPS(InitState(sites,"Up"));

  auto b = 4;

  SECTION("0-site")
    {
    auto Hpsi = LocalMPO(H,{"NumCenter=",0});

    CHECK(Hpsi.numCenter()==0);

    psi.position(b-1);

    auto l = linkIndex(psi,b-2);
    auto s = siteIndex(psi,b-1);
    auto [U,C] = polar(psi(b-1),{l,s},{"AddTags=","left","Prime=",0});
    psi.ref(b-1) = U;

    Hpsi.position(b,psi);

    auto phi = C;
    auto Hphi = phi;
    Hpsi.product(phi,Hphi);

    CHECK_CLOSE(norm(Hphi-noPrime(phi*Hpsi.L()*Hpsi.R())),0.);
    }

  SECTION("1-site")
    {
    auto Hpsi = LocalMPO(H,{"NumCenter=",1});

    CHECK(Hpsi.numCenter()==1);

    psi.position(b);
    Hpsi.position(b,psi);

    auto phi = psi(b);
    auto Hphi = phi;
    Hpsi.product(phi,Hphi);

    CHECK_CLOSE(norm(Hphi-noPrime(phi*Hpsi.L()*H(b)*Hpsi.R())),0.);
    }

  SECTION("2-site")
    {
    auto Hpsi = LocalMPO(H,{"NumCenter=",2});

    CHECK(Hpsi.numCenter()==2);

    psi.position(b);
    Hpsi.position(b,psi);

    auto phi = psi(b)*psi(b+1);
    auto Hphi = phi;
    Hpsi.product(phi,Hphi);

    CHECK_CLOSE(norm(Hphi-noPrime(phi*Hpsi.L()*H(b)*H(b+1)*Hpsi.R())),0.);
    }

  }

SECTION("LocalMPOSet")
  {
  int N = 10;
  auto sites = SpinHalf(N);
  auto ampo1 = AutoMPO(sites);
  auto ampo2 = AutoMPO(sites);
  for(int j = 1; j < N; ++j)
      {
      ampo1 += 0.5,"S+",j,"S-",j+1;
      ampo1 += 0.5,"S-",j,"S+",j+1;
      ampo2 +=     "Sz",j,"Sz",j+1;
      }
  auto H1 = toMPO(ampo1);
  auto H2 = toMPO(ampo2);
  auto H = std::vector<MPO>({H1,H2});
  auto psi = MPS(InitState(sites,"Up"));

  auto b = 4;

  SECTION("0-site")
    {
    auto Hpsi = LocalMPOSet(H,{"NumCenter=",0});

    CHECK(Hpsi.numCenter()==0);

    psi.position(b-1);

    auto l = linkIndex(psi,b-2);
    auto s = siteIndex(psi,b-1);
    auto [U,C] = polar(psi(b-1),{l,s},{"AddTags=","left","Prime=",0});
    psi.ref(b-1) = U;

    Hpsi.position(b,psi);

    auto phi = C;
    auto Hphi = phi;

    Hpsi.product(phi,Hphi);

    CHECK_CLOSE(norm(Hphi-noPrime(phi*Hpsi.L()[0]*Hpsi.R()[0]+
                                  phi*Hpsi.L()[1]*Hpsi.R()[1])),0.);
    }

  SECTION("1-site")
    {
    auto Hpsi = LocalMPOSet(H,{"NumCenter=",1});

    CHECK(Hpsi.numCenter()==1);

    psi.position(b);
    Hpsi.position(b,psi);

    auto phi = psi(b);
    auto Hphi = phi;
    Hpsi.product(phi,Hphi);

    CHECK_CLOSE(norm(Hphi-noPrime(phi*Hpsi.L()[0]*H[0](b)*Hpsi.R()[0]+
                                  phi*Hpsi.L()[1]*H[1](b)*Hpsi.R()[1])),0.);
    }

  SECTION("2-site")
    {
    auto Hpsi = LocalMPOSet(H,{"NumCenter=",2});

    CHECK(Hpsi.numCenter()==2);

    psi.position(b);
    Hpsi.position(b,psi);

    auto phi = psi(b)*psi(b+1);
    auto Hphi = phi;
    Hpsi.product(phi,Hphi);

    CHECK_CLOSE(norm(Hphi-noPrime(phi*Hpsi.L()[0]*H[0](b)*H[0](b+1)*Hpsi.R()[0]+
                                  phi*Hpsi.L()[1]*H[1](b)*H[1](b+1)*Hpsi.R()[1])),0.);
    }
  }
}


