//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOL_H
#define __ITENSOR_TEVOL_H

#include "mpo.h"
#include "bondgate.h"
#include "TEvolObserver.h"

namespace itensor {



//
// Evolves an MPS in real or imaginary time by an amount ttotal in steps
// of tstep using the list of bond gates provided.
//
// Arguments recognized:
//    "Verbose": if true, print useful information to stdout
//
template <class Iterable, class Tensor>
Real
gateTEvol(const Iterable& gatelist, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          const Args& args = Global::args());

template <class Iterable, class Tensor>
Real
gateTEvol(const Iterable& gatelist, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          Observer& obs,
          Args args = Global::args());

//
// Imaginary time evolve an MPS by an amount ttotal in time
// steps of tstep using the Hamiltonian MPO H.
//
// Arguments recognized:
//     "Verbose": print useful information to stdout
//     "Order": order at which to stop applying powers of H, 
//              setting order to p yields error of tstep^p
//     "Maxm": Maximum states kept each step
//     "Cutoff": Maximum truncation error each step
//     "Nsweep": Number of sweeps used to apply H to MPS (see fitApplyMPO)
//
template <class Tensor>
void
imagTEvol(const MPOt<Tensor>& H, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          const Args& args = Global::args());


//
//
// Implementations
//

template <class Iterable, class Tensor>
Real
gateTEvol(const Iterable& gatelist, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          Observer& obs,
          Args args)
    {
    const bool verbose = args.getBool("Verbose",false);
    const bool normalize = args.getBool("Normalize",true);

    const int nt = int(ttotal/tstep+(1e-9*(ttotal/tstep)));
    if(fabs(nt*tstep-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }

    Real tsofar = 0;
    Real tot_norm = psi.normalize();
    psi.position(gatelist.front().i());
    if(verbose) 
        {
        printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
        }
    for(int tt = 1; tt <= nt; ++tt)
        {
        for(const BondGate<Tensor>& G : gatelist)
            {
            const int lastpos = psi.orthoCenter();
            const int closest = abs(lastpos-G.i1()) < abs(lastpos-G.i2()) ? G.i1() : G.i2();
            psi.position(closest);
            applyGate(G,psi);
            }

        if(normalize)
            {
            tot_norm *= psi.normalize();
            }

        tsofar += tstep;

        args.add("TimeStepNum",tt);
        args.add("Time",tsofar);
        args.add("TotalTime",ttotal);
        obs.measure(args);
        }
    if(verbose) 
        {
        printfln("\nTotal time evolved = %.5f\n",tsofar);
        }

    return tot_norm;

    } // gateTEvol

template <class Iterable, class Tensor>
Real
gateTEvol(const Iterable& gatelist, 
          Real ttotal, 
          Real tstep, 
          MPSt<Tensor>& psi, 
          const Args& args)
    {
    TEvolObserver obs(args);
    return gateTEvol(gatelist,ttotal,tstep,psi,obs,args);
    }

} //namespace itensor


#endif
