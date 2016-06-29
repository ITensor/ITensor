#include "test.h"
#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

using namespace itensor;
using namespace std;

TEST_CASE("AutoMPO Test")
    {

SECTION("Longitudinal Ising")
    {
    int N = 10;
    SpinHalf sites(N);
    Real h = 0.2;

    AutoMPO ampo(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += "Sz",j,"Sz",j+1;
        }
    for(int j = 1; j <= N; ++j)
        {
        ampo += h,"Sz",j;
        }
    MPO H = toMPO<ITensor>(ampo);


    MPS psi(sites);
    for(int n = 1; n <= N; ++n) psi.Anc(n).randomize();

    Real E1 = psiHphi(psi,H,psi);

    Real E2 = 0.;
    for(int n = 1; n < N; ++n)
        {
        psi.position(n);
        ITensor op = sites.op("Sz",n)*sites.op("Sz",n+1);
        ITensor wf = psi.A(n)*psi.A(n+1);
        E2 += (dag(prime(wf,Site)) * op * wf).toReal();
        }
    for(int n = 1; n <= N; ++n)
        {
        psi.position(n);
        ITensor op = h*sites.op("Sz",n);
        ITensor wf = psi.A(n);
        E2 += (dag(prime(wf,Site)) * op * wf).toReal();
        }

    CHECK(fabs(E1-E2) < 1E-12);
    }

SECTION("Heisenberg")
    {
    int N = 10;
    SpinHalf sites(N);

    AutoMPO ampo(sites);
    for(int j = 1; j < N; ++j)
        {
        ampo += "Sz",j,"Sz",j+1;
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        }
    MPO H = toMPO<ITensor>(ampo);

    MPS psi(sites);
    for(int n = 1; n <= N; ++n) psi.Anc(n).randomize();

    Real E1 = psiHphi(psi,H,psi);

    Real E2 = 0.;
    for(int n = 1; n < N; ++n)
        {
        psi.position(n);
        ITensor op = sites.op("Sz",n)*sites.op("Sz",n+1);
        op += 0.5*sites.op("S+",n)*sites.op("S-",n+1);
        op += 0.5*sites.op("S-",n)*sites.op("S+",n+1);
        ITensor wf = psi.A(n)*psi.A(n+1);
        E2 += (dag(prime(wf,Site)) * op * wf).toReal();
        }

    CHECK(fabs(E1-E2) < 1E-12);
    }

    }
