#ifndef __HEISENBERG_H
#define __HEISENBERG_H

#include "itensor/mps/sites/spinhalf.h"
#include "itensor/mps/autompo.h"

namespace itensor {

MPO inline
Heisenberg(SpinHalf const& sites)
    {
    auto ampo = AutoMPO(sites);
    auto N = sites.N();
    for(int j = 1; j < N; ++j)
        {
        ampo += 0.5,"S+",j,"S-",j+1;
        ampo += 0.5,"S-",j,"S+",j+1;
        ampo +=     "Sz",j,"Sz",j+1;
        }
    auto H = MPO(ampo);

    return H;
    }

} //namespace itensor

#endif
