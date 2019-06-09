#include "itensor/all.h"

using namespace itensor;

int main()
    {
    auto Nx = 8,
         Ny = 4;
    auto N = Nx*Ny;
    auto sites = Electron(N);

    auto t = 1.0;
    auto U = 8.0;

    auto ampo = AutoMPO(sites);
    auto lattice = squareLattice(Nx,Ny,{"YPeriodic=",true});
    for(auto j : lattice)
        {
        ampo += -t,"Cdagup",j.s1,"Cup",j.s2;
        ampo += -t,"Cdagup",j.s2,"Cup",j.s1;
        ampo += -t,"Cdagdn",j.s1,"Cdn",j.s2;
        ampo += -t,"Cdagdn",j.s2,"Cdn",j.s1;
        }
    for(auto j : range1(N))
        {
        ampo += U,"Nupdn",j;
        }
    auto H = toMPO(ampo);

    auto state = InitState(sites);
    for(auto j : range1(N))
        {
        state.set(j,(j%2==1 ? "Up" : "Dn"));
        }
    auto psi0 = MPS(state);

    auto sweeps = Sweeps(5);
    sweeps.maxdim() = 10,20,100,200;
    sweeps.noise() = 1E-7,1E-8,1E-10,0;
    sweeps.cutoff() = 1E-6;
    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet=",true});

    return 0;
    }
