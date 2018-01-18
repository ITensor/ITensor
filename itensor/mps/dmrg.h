//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRG_H
#define __ITENSOR_DMRG_H

#include "itensor/eigensolver.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/localmpo_mps.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"


namespace itensor {

//
// Available DMRG methods:
//

//
//DMRG with an MPO
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Sweeps& sweeps,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and custom DMRGObserver
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Tensor& LH, const Tensor& RH,
     const Sweeps& sweeps,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
//and a custom observer
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Tensor& LH, const Tensor& RH,
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs,
     const Args& args = Global::args())
    {
    LocalMPO<Tensor> PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with a set of MPOs (lazily summed)
//(H vector is 0-indexed)
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const std::vector<MPOt<Tensor> >& Hset, 
     const Sweeps& sweeps,
     const Args& args = Global::args())
    {
    LocalMPOSet<Tensor> PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with a set of MPOs and a custom DMRGObserver
//(H vector is 0-indexed)
//
template <class Tensor>
Real 
dmrg(MPSt<Tensor>& psi, 
     const std::vector<MPOt<Tensor> >& Hset, 
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs,
     const Args& args = Global::args())
    {
    LocalMPOSet<Tensor> PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with a single Hamiltonian MPO and a set of 
//MPS to orthogonalize against
//(psis vector is 0-indexed)
//Named Args recognized:
// Weight - real number w > 0; calling dmrg(psi,H,psis,sweeps,Args("Weight",w))
//          sets the effective Hamiltonian to be
//          H + w * (|0><0| + |1><1| + ...) where |0> = psis[0], |1> = psis[1]
//          etc.
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const std::vector<MPSt<Tensor> >& psis, 
     const Sweeps& sweeps, 
     const Args& args = Global::args())
    {
    LocalMPO_MPS<Tensor> PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with a single Hamiltonian MPO, 
//a set of MPS to orthogonalize against, 
//and a custom DMRGObserver.
//(psis vector is 0-indexed)
//Named Args recognized:
// Weight - real number w > 0; calling dmrg(psi,H,psis,sweeps,Args("Weight",w))
//          sets the effective Hamiltonian to be
//          H + w * (|0><0| + |1><1| + ...) where |0> = psis[0], |1> = psis[1]
//          etc.
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const std::vector<MPSt<Tensor> >& psis, 
     const Sweeps& sweeps, 
     DMRGObserver<Tensor>& obs, 
     const Args& args = Global::args())
    {
    LocalMPO_MPS<Tensor> PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }



//
// DMRGWorker
//

template <class Tensor, class LocalOpT>
Real inline
DMRGWorker(MPSt<Tensor>& psi,
           LocalOpT& PH,
           const Sweeps& sweeps,
           const Args& args = Global::args())
    {
    DMRGObserver<Tensor> obs(psi,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

template <class Tensor, class LocalOpT>
Real
DMRGWorker(MPSt<Tensor>& psi,
           LocalOpT& PH,
           const Sweeps& sweeps,
           DMRGObserver<Tensor>& obs,
           Args args = Global::args())
    {
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? 0 : 1));

    const int N = psi.N();
    Real energy = NAN;

    psi.position(1);

    args.add("DebugLevel",debug_level);
    args.add("DoNormalize",true);
    
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        cpu_time sw_time;
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("Minm",sweeps.minm(sw));
        args.add("Maxm",sweeps.maxm(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

        if(!PH.doWrite()
           && args.defined("WriteM")
           && sweeps.maxm(sw) >= args.getInt("WriteM"))
            {
            if(!quiet)
                {
                println("\nTurning on write to disk, write_dir = ",
                        args.getString("WriteDir","./"));
                }

            //psi.doWrite(true);
            PH.doWrite(true);
            }

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(!quiet)
                {
                printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));
                }

            PH.position(b,psi);

            auto phi = psi.A(b)*psi.A(b+1);

            energy = davidson(PH,phi,args);
            
            auto spec = psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,args);


            if(!quiet)
                { 
                printfln("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d",
                          sweeps.cutoff(sw),
                          sweeps.minm(sw), 
                          sweeps.maxm(sw) );
                printfln("    Trunc. err=%.1E, States kept: %s",
                         spec.truncerr(),
                         showm(linkInd(psi,b)) );
                }

            obs.lastSpectrum(spec);

            args.add("AtBond",b);
            args.add("HalfSweep",ha);
            args.add("Energy",energy); 
            args.add("Truncerr",spec.truncerr()); 

            obs.measure(args);

            } //for loop over b

        auto sm = sw_time.sincemark();
        printfln("    Sweep %d CPU time = %s (Wall time = %s)",
                  sw,showtime(sm.time),showtime(sm.wall));

        if(obs.checkDone(args)) break;
    
        } //for loop over sw
    
    psi.normalize();

    return energy;
    }

} //namespace itensor


#endif
