#ifndef __J1J2_H
#define __J1J2_H

#include "itensor/all.h"

namespace itensor {

MPS inline
computeGroundState(SpinHalf const& sites, 
                   Real J2)
    {
    auto ampo = AutoMPO(sites);
    auto N = sites.N();
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    for(int j = 1; j < N-1; ++j)
        {
        ampo += 0.5*J2,"S+",j,"S-",j+2;
        ampo += 0.5*J2,"S-",j,"S+",j+2;
        ampo +=     J2,"Sz",j,"Sz",j+2;
        }
    auto H = MPO(ampo);

    auto psi = MPS(sites);

    auto sweeps = Sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-9;

    println("Starting ground state calculation for J2 = ",J2);

    dmrg(psi,H,sweeps,"Quiet");

    println("Done with ground state calculation.");

    return psi;
    }

} //namespace itensor

#endif
