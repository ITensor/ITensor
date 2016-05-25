#include "test.h"
#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/hubbard.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/util/print_macro.h"

using namespace itensor;

TEST_CASE("AutoMPO Test")
{

SECTION("Hubbard")
    {
    auto N = 10;
    auto sites = Hubbard(N);

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

    auto H = IQMPO(ampo);

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
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
        L = Vac;
        R = Vac;
        L.set(n,"Dn");
        R.set(n+1,"Dn");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
        }
    for(int n = 1; n < N; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n+1,"Up");
        R.set(n,"Up");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
        L = Vac;
        R = Vac;
        L.set(n+1,"Dn");
        R.set(n,"Dn");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
        }
    //
    // Check periodic t1
    //
    L = Vac;
    R = Vac;
    L.set(1,"Up");
    R.set(N,"Up");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
    L = Vac;
    R = Vac;
    L.set(N,"Up");
    R.set(1,"Up");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
    L = Vac;
    R = Vac;
    L.set(1,"Dn");
    R.set(N,"Dn");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);
    L = Vac;
    R = Vac;
    L.set(N,"Dn");
    R.set(1,"Dn");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t1);

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
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),+t1);

    //
    // Check t2
    //
    for(int n = 1; n < N-1; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n,"Up");
        R.set(n+2,"Up");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t2);
        }
    for(int n = 1; n < N-1; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n+2,"Up");
        R.set(n,"Up");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t2);
        }
    //
    // Check periodic t2
    //
    L = Vac;
    R = Vac;
    L.set(2,"Up");
    R.set(N,"Up");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t2);
    L = Vac;
    R = Vac;
    L.set(N,"Up");
    R.set(2,"Up");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t2);
    L = Vac;
    R = Vac;
    L.set(2,"Dn");
    R.set(N,"Dn");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t2);
    L = Vac;
    R = Vac;
    L.set(N,"Dn");
    R.set(2,"Dn");
    CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),-t2);

    //
    // Check U
    //
    for(int n = 1; n <= N; ++n)
        {
        L = Vac;
        R = Vac;
        L.set(n,"UpDn");
        R.set(n,"UpDn");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),U);
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
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),V1);
        L = Vac;
        R = Vac;
        L.set(n,"Up");
        L.set(n+1,"Dn");
        R.set(n,"Up");
        R.set(n+1,"Dn");
        CHECK_CLOSE(overlap(IQMPS(L),H,IQMPS(R)),V1);
        }

    }

SECTION("No QN MPO")
    {
    auto N = 10;
    auto h = 0.732;
    auto sites = SpinHalf(N);
    auto ampo = AutoMPO(sites);
    for(auto j = 1; j < N; ++j)
        {
        ampo += "Sx",j,"Sx",j+1;
        }
    for(auto j = 1; j <= N; ++j)
        {
        ampo += h,"Sz",j;
        }
    auto H = MPO(ampo);

    auto AllUp = InitState(sites,"Up");
    auto L = AllUp;
    auto R = AllUp;
    CHECK_CLOSE(overlap(MPS(L),H,MPS(R)),N*h/2.);

    L = AllUp;
    R = AllUp;
    L.set(1,"Dn");
    R.set(1,"Dn");
    CHECK_CLOSE(overlap(MPS(L),H,MPS(R)),(N-2)*h/2.);

    L = AllUp;
    R = AllUp;
    L.set(1,"Dn");
    L.set(2,"Dn");
    CHECK_CLOSE(overlap(MPS(L),H,MPS(R)),1./4.);

    L = AllUp;
    R = AllUp;
    L.set(3,"Dn");
    R.set(4,"Dn");
    CHECK_CLOSE(overlap(MPS(L),H,MPS(R)),1./4.);
    }

}
