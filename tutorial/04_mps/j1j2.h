#ifndef __J1J2_H
#define __J1J2_H

#include "core.h"
#include "sites/spinhalf.h"
#include "hams/J1J2Chain.h"

using namespace itensor;

MPS inline
computeGroundState(const SpinHalf& model, Real J2)
    {
    MPO H = J1J2Chain(model,Opt("J2",J2));

    MPS psi(model);

    Sweeps sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-9;

    println("Starting ground state calculation for J2 = ",J2);

    dmrg(psi,H,sweeps,"Quiet");

    println("Done with ground state calculation.");

    return psi;
    }

#endif
