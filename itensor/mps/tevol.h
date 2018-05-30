//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEVOL_H
#define __ITENSOR_TEVOL_H

#include "itensor/mps/mpo.h"
#include "itensor/mps/bondgate.h"
#include "itensor/mps/TEvolObserver.h"

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
//
// Implementations
//

template <class Iterable, class Tensor>
Real
gateTEvol(Iterable const& gatelist, 
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

    if(verbose) 
        {
        printfln("Taking %d steps of timestep %.5f, total time %.5f",nt,tstep,ttotal);
        }

    psi.position(gatelist.front().i1());
    Real tot_norm = norm(psi);

    Real tsofar = 0;
    for(int tt = 1; tt <= nt; ++tt)
        {
        auto g = gatelist.begin();
        while(g != gatelist.end())
            {
            auto i1 = g->i1();
            auto i2 = g->i2();
            auto AA = psi.A(i1)*psi.A(i2)*g->gate();
            AA.mapprime(1,0,Site);

            ++g;
            if(g != gatelist.end())
                {
                //Look ahead to next gate position
                auto ni1 = g->i1();
                auto ni2 = g->i2();
                //SVD AA to restore MPS form
                //before applying current gate
                if(ni1 >= i2)
                    {
                    psi.svdBond(i1,AA,Fromleft,args);
                    psi.position(ni1); //does no work if position already ni1
                    }
                else
                    {
                    psi.svdBond(i1,AA,Fromright,args);
                    psi.position(ni2); //does no work if position already ni2
                    }
                }
            else
                {
                //No next gate to analyze, just restore MPS form
                psi.svdBond(i1,AA,Fromright,args);
                }
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
