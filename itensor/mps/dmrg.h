//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRG_H
#define __ITENSOR_DMRG_H

#include "itensor/iterativesolvers.h"
#include "itensor/mps/localmposet.h"
#include "itensor/mps/localmpo_mps.h"
#include "itensor/mps/sweeps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/cputime.h"


namespace itensor {

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           Args const& args = Args::global());

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           DMRGObserver & obs,
           Args args = Args::global());

//
// Available DMRG methods:
//

//
//DMRG with an MPO
//
Real inline
dmrg(MPS & psi, 
     MPO const& H, 
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPO PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and custom DMRGObserver
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     Sweeps const& sweeps, 
     DMRGObserver & obs,
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPO PH(H,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
// LH - H1 - H2 - ... - HN - RH
//(ok if one or both of LH, RH default constructed)
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     ITensor const& LH, 
     ITensor const& RH,
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPO PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with an MPO and boundary tensors LH, RH
//and a custom observer
//
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     ITensor const& LH, 
     ITensor const& RH,
     Sweeps const& sweeps, 
     DMRGObserver& obs,
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPO PH(H,LH,RH,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

//
//DMRG with a set of MPOs (lazily summed)
//(H vector is 0-indexed)
//
Real inline
dmrg(MPS& psi, 
     std::vector<MPO> const& Hset, 
     Sweeps const& sweeps,
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPOSet PH(Hset,args);
    Real energy = DMRGWorker(psi,PH,sweeps,args);
    return energy;
    }

//
//DMRG with a set of MPOs and a custom DMRGObserver
//(H vector is 0-indexed)
//
Real inline
dmrg(MPS& psi, 
     std::vector<MPO> const& Hset, 
     Sweeps const& sweeps, 
     DMRGObserver& obs,
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPOSet PH(Hset,args);
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
Real inline
dmrg(MPS& psi, 
     MPO const& H, 
     std::vector<MPS> const& psis, 
     Sweeps const& sweeps, 
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPO_MPS PH(H,psis,args);
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
Real inline
dmrg(MPS & psi, 
     MPO const& H, 
     std::vector<MPS> const& psis, 
     Sweeps const& sweeps, 
     DMRGObserver& obs, 
     Args const& args = Args::global())
    {
    if(args.defined("Maxm"))
      Error("Error in dmrg: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in dmrg: Arg Minm is deprecated in favor of MinDim.");

    LocalMPO_MPS PH(H,psis,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }



//
// DMRGWorker
//

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           Args const& args)
    {
    if(args.defined("Maxm"))
      Error("Error in DMRGWorker: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in DMRGWorker: Arg Minm is deprecated in favor of MinDim.");

    DMRGObserver obs(psi,args);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,args);
    return energy;
    }

template<class LocalOpT>
Real
DMRGWorker(MPS & psi,
           LocalOpT & PH,
           Sweeps const& sweeps,
           DMRGObserver & obs,
           Args args)
    {
    if(args.defined("Maxm"))
      Error("Error in DMRGWorker: Arg Maxm is deprecated in favor of MaxDim.");
    if(args.defined("Minm"))
      Error("Error in DMRGWorker: Arg Minm is deprecated in favor of MinDim.");
 
    const bool quiet = args.getBool("Quiet",false);
    const int debug_level = args.getInt("DebugLevel",(quiet ? 0 : 1));

    const int N = length(psi);
    Real energy = NAN;

    psi.position(1);

    args.add("DebugLevel",debug_level);
    args.add("DoNormalize",true);
    
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        cpu_time sw_time;
        args.add("Sweep",sw);
        args.add("NSweep",sweeps.nsweep());
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

        if(!PH.doWrite()
           && args.defined("WriteM")
           && sweeps.maxdim(sw) >= args.getInt("WriteM"))
            {
            if(!quiet)
                {
                println("\nTurning on write to disk, write_dir = ",
                        args.getString("WriteDir","./"));
                }

            //psi.doWrite(true);
            PH.doWrite(true,args);
            }

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(!quiet)
                {
                printfln("Sweep=%d, HS=%d, Bond=%d/%d",sw,ha,b,(N-1));
                }

            PH.position(b,psi);

            auto phi = psi(b)*psi(b+1);

            energy = davidson(PH,phi,args);
            
            auto spec = psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,args);


            if(!quiet)
                { 
                printfln("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d",
                          sweeps.cutoff(sw),
                          sweeps.mindim(sw), 
                          sweeps.maxdim(sw) );
                printfln("    Trunc. err=%.1E, States kept: %s",
                         spec.truncerr(),
                         showDim(linkInd(psi,b)) );
                }

            obs.lastSpectrum(spec);

            args.add("AtBond",b);
            args.add("HalfSweep",ha);
            args.add("Energy",energy); 
            args.add("Truncerr",spec.truncerr()); 

            obs.measure(args);

            } //for loop over b

        auto sm = sw_time.sincemark();
        printfln("    Sweep %d/%d CPU time = %s (Wall time = %s)",
                  sw,sweeps.nsweep(),showtime(sm.time),showtime(sm.wall));

        if(obs.checkDone(args)) break;
    
        } //for loop over sw
    
    psi.normalize();

    return energy;
    }

} //namespace itensor


#endif
