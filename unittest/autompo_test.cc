#include "test.h"
#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/electron.h"
#include "itensor/mps/sites/fermion.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/print_macro.h"

#include "ExpIsing.h"
#include "ExpHeisenberg.h"

using namespace itensor;

TEST_CASE("AutoMPO Test")
{

SECTION("Electron")
    {
    auto N = 10;
    auto sites = Electron(N);

    auto U = 0.1;
    auto t1 = 0.7;
    auto t2 = 0.25;
    auto V1 = 2.01;

    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= N; ++i)
        {
        ampo += U,"Nupdn",i;
        }
    for(int b = 1; b < N; ++b)
        {
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;
        ampo += V1,"Ntot",b,"Ntot",b+1;
        }
    ampo += -t1,"Cdagup",1,"Cup",N;
    ampo += -t1,"Cdagup",N,"Cup",1;
    ampo += -t1,"Cdagdn",1,"Cdn",N;
    ampo += -t1,"Cdagdn",N,"Cdn",1;

    for(int b = 1; b < N-1; ++b)
        {
        ampo += -t2,"Cdagup",b,"Cup",b+2;
        ampo += -t2,"Cdagup",b+2,"Cup",b;
        ampo += -t2,"Cdagdn",b,"Cdn",b+2;
        //Change order of this one to
        //make sure correct fermion sign
        //is put in:
        ampo += +t2,"Cdn",b,"Cdagdn",b+2;
        }
    ampo += -t2,"Cdagup",2,"Cup",N;
    ampo += -t2,"Cdagup",N,"Cup",2;
    ampo += -t2,"Cdagdn",2,"Cdn",N;
    ampo += -t2,"Cdagdn",N,"Cdn",2;

    ampo += -t2,"Cdagup",1,"Cup",N-1;
    ampo += -t2,"Cdagup",N-1,"Cup",1;
    ampo += -t2,"Cdagdn",1,"Cdn",N-1;
    ampo += -t2,"Cdagdn",N-1,"Cdn",1;

    auto H = toMPO(ampo);

    auto Vac = InitState(sites,"Emp");
    auto L = Vac;
    auto R = Vac;

    //
    // Check t1
    //
    for(int n = 1; n < N; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n,"Up");
        R.set(n+1,"Up");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
        L = Vac;
        R = Vac;
        L.set(n,"Dn");
        R.set(n+1,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
        }
    for(int n = 1; n < N; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n+1,"Up");
        R.set(n,"Up");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
        L = Vac;
        R = Vac;
        L.set(n+1,"Dn");
        R.set(n,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
        }
    //
    // Check periodic t1
    //
    L = Vac;
    R = Vac;
    L.set(1,"Up");
    R.set(N,"Up");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
    L = Vac;
    R = Vac;
    L.set(N,"Up");
    R.set(1,"Up");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
    L = Vac;
    R = Vac;
    L.set(1,"Dn");
    R.set(N,"Dn");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);
    L = Vac;
    R = Vac;
    L.set(N,"Dn");
    R.set(1,"Dn");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t1);

    //
    // Check that periodic t1 is fermionic
    //
    L = Vac;
    R = Vac;
    L.set(N,"Dn");
    R.set(1,"Dn");
    //Put particle in between that the other
    //has to "hop" over:
    L.set(2,"Up");
    R.set(2,"Up");
    //Should change the sign of resulting matrix element:
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),+t1);

    //
    // Check t2
    //
    for(int n = 1; n < N-1; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n,"Up");
        R.set(n+2,"Up");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t2);
        }
    for(int n = 1; n < N-1; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n+2,"Up");
        R.set(n,"Up");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t2);
        }
    //
    // Check periodic t2
    //
    L = Vac;
    R = Vac;
    L.set(2,"Up");
    R.set(N,"Up");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t2);
    L = Vac;
    R = Vac;
    L.set(N,"Up");
    R.set(2,"Up");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t2);
    L = Vac;
    R = Vac;
    L.set(2,"Dn");
    R.set(N,"Dn");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t2);
    L = Vac;
    R = Vac;
    L.set(N,"Dn");
    R.set(2,"Dn");
    CHECK_CLOSE(inner(MPS(L),H,MPS(R)),-t2);

    //
    // Check U
    //
    for(int n = 1; n <= N; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n,"UpDn");
        R.set(n,"UpDn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),U);
        }

    //
    // Check V1
    //
    for(int n = 1; n < N; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n,"Up");
        L.set(n+1,"Up");
        R.set(n,"Up");
        R.set(n+1,"Up");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),V1);
        L = Vac;
        R = Vac;
        L.set(n,"Up");
        L.set(n+1,"Dn");
        R.set(n,"Up");
        R.set(n+1,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),V1);
        }

    }

SECTION("No QN MPO")
    {
    auto N = 10;
    auto h = 0.732;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto ampo = AutoMPO(sites);
    for(auto j = 1; j < N; ++j)
        {
        ampo += "Sx",j,"Sx",j+1;
        }
    for(auto j = 1; j <= N; ++j)
        {
        ampo += h,"Sz",j;
        }

    SECTION("Exact version")
        {
        auto H = toMPO(ampo,{"Exact=",true});

        auto AllUp = InitState(sites,"Up");
        auto L = AllUp;
        auto R = AllUp;
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),N*h/2.);

        L = AllUp;
        R = AllUp;
        L.set(1,"Dn");
        R.set(1,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),(N-2)*h/2.);

        L = AllUp;
        R = AllUp;
        L.set(1,"Dn");
        L.set(2,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),1./4.);

        L = AllUp;
        R = AllUp;
        L.set(3,"Dn");
        R.set(4,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),1./4.);
        }

    SECTION("Approx version")
        {
        auto H = toMPO(ampo,{"Exact",false});

        auto AllUp = InitState(sites,"Up");
        auto L = AllUp;
        auto R = AllUp;
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),N*h/2.);

        L = AllUp;
        R = AllUp;
        L.set(1,"Dn");
        R.set(1,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),(N-2)*h/2.);

        L = AllUp;
        R = AllUp;
        L.set(1,"Dn");
        L.set(2,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),1./4.);

        L = AllUp;
        R = AllUp;
        L.set(3,"Dn");
        R.set(4,"Dn");
        CHECK_CLOSE(inner(MPS(L),H,MPS(R)),1./4.);
        }
    }

SECTION("Single Site Ops")
    {
    int L = 10;
    auto sites = Electron(L);
    auto ampo = AutoMPO(sites);

    SECTION("Diagonal op")
        {
        auto n = 4;
        ampo += "Nup",n;
        auto Op = toMPO(ampo);
        for(auto i : range1(L))
            {
            auto state = InitState(sites,"Emp");
            state.set(i,"Up");
            auto psi = MPS(state);
            auto x = inner(psi,Op,psi);
            if(i == n) CHECK_CLOSE(x,1.);
            else       CHECK_CLOSE(x,0.);
            }
        }

    //SECTION("Off diagonal, fermionic op")
    //    {
    //    auto n = 4;
    //    ampo += "Cdagup",n;
    //    auto Op = toMPO(ampo);
    //    //PrintData(Op(1));
    //    //PrintData(Op(2));
    //    //PrintData(Op(3));
    //    //PrintData(Op(4));
    //    for(auto i : range1(L))
    //        {
    //        auto state = InitState(sites,"Emp");
    //        auto vac = MPS(state);
    //        state.set(i,"Up");
    //        auto psi = MPS(state);
    //        auto x = inner(psi,Op,vac);
    //        if(i == n) CHECK_CLOSE(x,1.);
    //        else       CHECK_CLOSE(x,0.);
    //        }
    //    }
    }

SECTION("toExpH ITensor (no QNs)")
    {
    int N = 10;
    Real h = 0.2;
    Real tau = 0.01234;

    auto sites = SpinHalf(N,{"ConserveQNs=",false}); //make a chain of N spin 1/2's

    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += "Sz",j,"Sz",j+1;
        }
    for(int j = 1; j <= N; ++j)
        {
        ampo += -h,"Sx",j;
        }

    SECTION("Real time")
        {
        auto expH = toExpH(ampo,tau*1_i);
        auto expHexact = MPO(ExpIsing(sites,tau*1_i,{"h",h}));
        auto psi = randomMPS(sites);
        auto xpsi = applyMPO(expHexact,psi,{"Method","DensityMatrix"});
        auto xnrm2 = real(innerC(xpsi,xpsi));
        auto apsi = applyMPO(expH,psi,{"Method","DensityMatrix"});
        auto anrm2 = real(innerC(apsi,apsi));
        CHECK_CLOSE(innerC(xpsi,apsi)/sqrt(xnrm2*anrm2),1.);
      }
    SECTION("Imaginary time")
        {
        auto expH = toExpH(ampo,tau);
        auto expHexact = MPO(ExpIsing(sites,tau,{"h",h}));
        auto psi = randomMPS(sites);
        auto xpsi = applyMPO(expHexact,psi,{"Method","DensityMatrix"});
        auto xnrm2 = inner(xpsi,xpsi);
        auto apsi = applyMPO(expH,psi,{"Method","DensityMatrix"});
        auto anrm2 = inner(apsi,apsi);
        CHECK_CLOSE(inner(xpsi,apsi)/sqrt(xnrm2*anrm2),1.);
        }
  }

SECTION("toExpH ITensor (QN conservation)")
    {
    int N = 10;
    Real tau = 0.01;

    auto sites = SpinHalf(N); //make a chain of N spin 1/2's

    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N; ++j)
    {
    ampo +=     "Sz",j,"Sz",j+1;
    ampo += 0.5,"S+",j,"S-",j+1;
    ampo += 0.5,"S-",j,"S+",j+1;
    }

    auto state = InitState(sites);
    for(auto j : range1(N)) state.set(j,j%2==0 ? "Up" : "Dn");
    auto psi = MPS(state);

    SECTION("Real time")
        {
        auto expH = toExpH(ampo,tau*1_i);
        auto expHexact = MPO(ExpHeisenberg(sites,tau*1_i));
        auto xpsi = applyMPO(expHexact,psi,{"Method","DensityMatrix"});
        auto xnrm2 = real(innerC(xpsi,xpsi));
        auto apsi = applyMPO(expH,psi,{"Method","DensityMatrix"});
        auto anrm2 = real(innerC(apsi,apsi));
        CHECK_CLOSE(innerC(xpsi,apsi)/sqrt(xnrm2*anrm2),1.0);
        }
    SECTION("Imaginary time")
        {
        auto expH = toExpH(ampo,tau);
        auto expHexact = MPO(ExpHeisenberg(sites,tau));
        auto xpsi = applyMPO(expHexact,psi,{"Method","DensityMatrix"});
        auto xnrm2 = inner(xpsi,xpsi);
        auto apsi = applyMPO(expH,psi,{"Method","DensityMatrix"});
        auto anrm2 = inner(xpsi,xpsi);
        CHECK_CLOSE(inner(xpsi,apsi)/sqrt(xnrm2*anrm2),1.0);
        }
    }

SECTION("Electron, Complex Hopping")
    {
    auto N = 10;
    auto sites = Electron(N);

    auto t = 0.5;
    auto phi = 0.123;

    auto ampo = AutoMPO(sites);
    for(auto b : range1(N-1))
        {
        ampo += (-t*std::polar(1.,-phi/4.)),"Cdagup",b,"Cup",b+1;
        ampo += (-t*std::polar(1.,+phi/4.)),"Cdagup",b+1,"Cup",b;
        ampo += (-t*std::polar(1.,-phi/4.)),"Cdagdn",b,"Cdn",b+1;
        ampo += (-t*std::polar(1.,+phi/4.)),"Cdagdn",b+1,"Cdn",b;
        }
    auto H = toMPO(ampo);

    for(auto b : range1(N-1))
        {
        auto rstate1u = InitState(sites,"Emp");
        rstate1u.set(b,"Up");
        auto rstate2u = InitState(sites,"Emp");
        rstate2u.set(b+1,"Up");
        auto rpsi1u = MPS(rstate1u);
        auto rpsi2u = MPS(rstate2u);
        CHECK_CLOSE(innerC(rpsi2u,H,rpsi1u),-t*std::polar(1.,+phi/4.));

        auto rstate1d = InitState(sites,"Emp");
        rstate1d.set(b,"Dn");
        auto rstate2d = InitState(sites,"Emp");
        rstate2d.set(b+1,"Dn");
        auto rpsi1d = MPS(rstate1d);
        auto rpsi2d = MPS(rstate2d);
        CHECK_CLOSE(innerC(rpsi2d,H,rpsi1d),-t*std::polar(1.,+phi/4.));
        }

    for(auto c : range1(2,N))
        {
        auto lstate1u = InitState(sites,"Emp");
        lstate1u.set(c,"Up");
        auto lstate2u = InitState(sites,"Emp");
        lstate2u.set(c-1,"Up");
        auto lpsi1u = MPS(lstate1u);
        auto lpsi2u = MPS(lstate2u);
        CHECK_CLOSE(innerC(lpsi2u,H,lpsi1u),-t*std::polar(1.,-phi/4.));

        auto lstate1d = InitState(sites,"Emp");
        lstate1d.set(c,"Up");
        auto lstate2d = InitState(sites,"Emp");
        lstate2d.set(c-1,"Up");
        auto lpsi1d = MPS(lstate1d);
        auto lpsi2d = MPS(lstate2d);
        CHECK_CLOSE(innerC(lpsi2d,H,lpsi1d),-t*std::polar(1.,-phi/4.));
        }
    }

SECTION("Electron with no QNs")
    {
    int N = 2;
    auto t0 = 5.0;

    auto sites = Fermion(N, {"ConserveQNs=",false});

    auto ampo = AutoMPO(sites);

    ampo += -t0,"Cdag",1,"C",2;
    ampo += -t0,"Cdag",2,"C",1;

    auto H = toMPO(ampo);

    auto HT = H(1)*H(2);
    auto HTc = dag(swapPrime(HT,0,1));
    CHECK(norm(HT-HTc) < 1E-10);
    }

SECTION("Ladder with Complex Hopping")
    {
    auto N = 8;
    auto sites = Electron(N);

    auto tpara = 0.5;
    auto tperp = 0.2;
    auto phi = 0.123;

    auto ampo = AutoMPO(sites);
	for(int j=1; j<=N; j+=1)
		{
        if(j < N-1)
            {
            if(phi==0)
                {
                ampo += -tpara,"Cdagup",j,"Cup",j+2;
                ampo += -tpara,"Cdagup",j+2,"Cup",j;
                ampo += -tpara,"Cdagdn",j,"Cdn",j+2;
                ampo += -tpara,"Cdagdn",j+2,"Cdn",j;
                }
            else
                {
                ampo += (-(tpara)*std::polar(1.,-phi/2.)),"Cdagup",j,"Cup",j+2;
                ampo += (-(tpara)*std::polar(1.,+phi/2.)),"Cdagup",j+2,"Cup",j;
                ampo += (-(tpara)*std::polar(1.,-phi/2.)),"Cdagdn",j,"Cdn",j+2;
                ampo += (-(tpara)*std::polar(1.,+phi/2.)),"Cdagdn",j+2,"Cdn",j;
                }
            }
        if(j%2 == 1)
            {
            ampo += -tperp,"Cdagup",j,"Cup",j+1;
            ampo += -tperp,"Cdagup",j+1,"Cup",j;
            ampo += -tperp,"Cdagdn",j,"Cdn",j+1;
            ampo += -tperp,"Cdagdn",j+1,"Cdn",j;
            }
		}
	auto Hx = toMPO(ampo,{"Exact",true});
	auto Ha = toMPO(ampo,{"Exact",false});

	for(int j=1; j<=N; j+=1)
		{
        if(j < N-1)
            {
            //Right hop
            {
            auto s1 = InitState(sites,"Emp");
            s1.set(j,"Up");
            auto s2 = InitState(sites,"Emp");
            s2.set(j+2,"Up");
            auto p1 = MPS(s1);
            auto p2 = MPS(s2);
            CHECK_CLOSE(innerC(p2,Hx,p1),-(tpara)*std::polar(1.,+phi/2.));
            CHECK_CLOSE(innerC(p2,Ha,p1),-(tpara)*std::polar(1.,+phi/2.));
            }

            //Left hop
            {
            auto s1 = InitState(sites,"Emp");
            s1.set(j+2,"Up");
            auto s2 = InitState(sites,"Emp");
            s2.set(j,"Up");
            auto p1 = MPS(s1);
            auto p2 = MPS(s2);
            CHECK_CLOSE(innerC(p2,Hx,p1),-(tpara)*std::polar(1.,-phi/2.));
            CHECK_CLOSE(innerC(p2,Ha,p1),-(tpara)*std::polar(1.,-phi/2.));
            }

            }
        if(j%2 == 1)
            {
            //Right hop
            {
            auto s1 = InitState(sites,"Emp");
            s1.set(j,"Up");
            auto s2 = InitState(sites,"Emp");
            s2.set(j+1,"Up");
            auto p1 = MPS(s1);
            auto p2 = MPS(s2);
            CHECK_CLOSE(innerC(p2,Hx,p1),-tperp);
            CHECK_CLOSE(innerC(p2,Ha,p1),-tperp);
            }

            //Left hop
            {
            auto s1 = InitState(sites,"Emp");
            s1.set(j+1,"Up");
            auto s2 = InitState(sites,"Emp");
            s2.set(j,"Up");
            auto p1 = MPS(s1);
            auto p2 = MPS(s2);
            CHECK_CLOSE(innerC(p2,Hx,p1),-tperp);
            CHECK_CLOSE(innerC(p2,Ha,p1),-tperp);
            }
            }
        }
    }

SECTION("Fermion")
    {
    auto N = 4;
    auto sites = Fermion(N);
    auto ampo = AutoMPO(sites);
    auto t = 0.5;
    for(auto b : range1(N-1))
        {
        ampo += -t,"Cdag",b,"C",b+1;
        ampo += -t,"Cdag",b+1,"C",b;
        }
    //Approx MPO construction
    auto Ha = toMPO(ampo,{"Exact",false});
    //Exact MPO construction
    auto Hx = toMPO(ampo,{"Exact",true});

    for(auto j : range1(N))
        {
        CHECK(not isComplex(Ha(j)));
        }

    for(auto b : range1(N-1))
        {
        auto rstate1 = InitState(sites,"Emp");
        rstate1.set(b,"Occ");
        auto rstate2 = InitState(sites,"Emp");
        rstate2.set(b+1,"Occ");
        auto rpsi1 = MPS(rstate1);
        auto rpsi2 = MPS(rstate2);
        CHECK_CLOSE(inner(rpsi2,Hx,rpsi1),-t);
        CHECK_CLOSE(inner(rpsi2,Ha,rpsi1),-t);
        }

    for(auto b : range1(2,N))
        {
        auto lstate1 = InitState(sites,"Emp");
        lstate1.set(b,"Occ");
        auto lstate2 = InitState(sites,"Emp");
        lstate2.set(b-1,"Occ");
        auto lpsi1 = MPS(lstate1);
        auto lpsi2 = MPS(lstate2);
        CHECK_CLOSE(inner(lpsi2,Hx,lpsi1),-t);
        CHECK_CLOSE(inner(lpsi2,Ha,lpsi1),-t);
        }
    }

SECTION("Mixed Fermion and Non-Fermion Sites")
    {
    // This test checks whether fermionic and non-fermionic
    // sites can be successfully mixed together and used
    // in the AutoMPO system. One possible issue is whether
    // the Jordan-Wigner "F" operator is defined for
    // non-fermionic sites
    using FSpin = MixedSiteSet<FermionSite,SpinHalfSite>;
    int N = 6;
    auto t = 1.0;

    auto sites = FSpin(N);

    auto ampo = AutoMPO(sites);
    for(int j = 1; j < N-1; j += 2)
        {
        ampo += -t,"Cdag",j,"C",j+2;
        ampo += -t,"Cdag",j+2,"C",j;
        }

    MPO H;
    CHECK_NOTHROW(H = toMPO(ampo));
    }

}
