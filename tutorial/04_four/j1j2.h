#ifndef __J1J2_H
#define __J1J2_H

#include "core.h"
#include "model/spinhalf.h"
#include "hams/J1J2Chain.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format


MPS inline
computeGroundState(const SpinHalf& model, Real J2)
    {
    MPO H = J1J2Chain(model,Opt("J2",J2));

    MPS psi(model);

    Sweeps sweeps(5);
    sweeps.maxm() = 50,50,100,100,200;
    sweeps.cutoff() = 1E-9;

    Cout << Format("Starting ground state calculation for J2 = %.3f")  
            % J2 << Endl;

    dmrg(psi,H,sweeps,Opt("Quiet"));

    Cout << "Done with ground state calculation." << Endl;

    return psi;
    }



#undef Cout
#undef Endl
#undef Format

#endif
