//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRG_H
#define __ITENSOR_DMRG_H

#include "eigensolver.h"
#include "localmposet.h"
#include "localmpo_mps.h"
#include "Sweeps.h"
#include "DMRGObserver.h"

#define Cout std::cout
#define Endl std::endl
#define Format boost::format

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
     const OptSet& opts = Global::opts())
    {
    LocalMPO<Tensor> PH(H,opts);
    Real energy = DMRGWorker(psi,PH,sweeps,opts);
    return energy;
    }

//
//DMRG with an MPO and custom Observer
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const Sweeps& sweeps, 
     Observer& obs,
     const OptSet& opts = Global::opts())
    {
    LocalMPO<Tensor> PH(H,opts);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,opts);
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
     const OptSet& opts = Global::opts())
    {
    LocalMPOSet<Tensor> PH(Hset,opts);
    Real energy = DMRGWorker(psi,PH,sweeps,opts);
    return energy;
    }

//
//DMRG with a set of MPOs and a custom Observer
//(H vector is 0-indexed)
//
template <class Tensor>
Real 
dmrg(MPSt<Tensor>& psi, 
     const std::vector<MPOt<Tensor> >& Hset, 
     const Sweeps& sweeps, 
     Observer& obs,
     const OptSet& opts = Global::opts())
    {
    LocalMPOSet<Tensor> PH(Hset,opts);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,opts);
    return energy;
    }

//
//DMRG with a single Hamiltonian MPO and a set of 
//MPS to orthogonalize against
//(psis vector is 0-indexed)
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const std::vector<MPSt<Tensor> >& psis, 
     const Sweeps& sweeps, 
     const OptSet& opts = Global::opts())
    {
    LocalMPO_MPS<Tensor> PH(H,psis,opts);
    Real energy = DMRGWorker(psi,PH,sweeps,opts);
    return energy;
    }

//
//DMRG with a single Hamiltonian MPO, 
//a set of MPS to orthogonalize against, 
//and a custom Observer.
//(psis vector is 0-indexed)
//
template <class Tensor>
Real
dmrg(MPSt<Tensor>& psi, 
     const MPOt<Tensor>& H, 
     const std::vector<MPSt<Tensor> >& psis, 
     const Sweeps& sweeps, Observer& obs, 
     const OptSet& opts = Global::opts())
    {
    LocalMPO_MPS<Tensor> PH(H,psis,opts);
    Real energy = DMRGWorker(psi,PH,sweeps,obs,opts);
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
           const OptSet& opts = Global::opts())
    {
    //Default Observer type is DMRGObserver
    DMRGObserver obs;
    Real energy = DMRGWorker(psi,PH,sweeps,obs,opts);
    return energy;
    }

template <class Tensor, class LocalOpT>
Real
DMRGWorker(MPSt<Tensor>& psi,
           LocalOpT& PH,
           const Sweeps& sweeps,
           Observer& obs,
           const OptSet& opts = Global::opts())
    {
    const Real orig_cutoff = psi.cutoff(),
               orig_noise  = psi.noise();
    const int orig_minm = psi.minm(), 
              orig_maxm = psi.maxm();

    const bool quiet = opts.getBool("Quiet",false);
    const int debug_level = opts.getInt("DebugLevel",(quiet ? 0 : 1));

    const int N = psi.N();
    Real energy = NAN;

    psi.position(1);
    
    Eigensolver solver(opts & Opt("DebugLevel",debug_level));

    const Opt doNorm = DoNormalize(true);
    
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        psi.cutoff(sweeps.cutoff(sw)); 
        psi.minm(sweeps.minm(sw)); 
        psi.maxm(sweeps.maxm(sw));
        psi.noise(sweeps.noise(sw));
        solver.maxIter(sweeps.niter(sw));

        if(!PH.doWrite() &&
            Global::opts().defined("WriteM") &&
            sweeps.maxm(sw) >= Global::opts().getInt("WriteM"))
            {
            std::string write_dir = Global::opts().getString("WriteDir","./");

            if(!quiet)
                {
                Cout << "\nTurning on write to disk, write_dir = " 
                     << write_dir << Endl;
                }

            psi.doWrite(true);
            PH.doWrite(true);
            }

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(!quiet)
                {
                Cout << Format("Sweep=%d, HS=%d, Bond=(%d,%d)") 
                        % sw % ha % b % (b+1) << Endl;
                }

            PH.position(b,psi);

            Tensor phi = psi.bondTensor(b);

            energy = solver.davidson(PH,phi);
            
            psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,doNorm);

            if(!quiet)
                { 
                Cout << 
                    Format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d") 
                        % sweeps.cutoff(sw) 
                        % sweeps.minm(sw) 
                        % sweeps.maxm(sw) 
                        << Endl;
                Cout << Format("    Trunc. err=%.1E, States kept=%s")
                        % psi.svd().truncerr(b) 
                        % showm(psi.LinkInd(b)) 
                        << Endl;
                }

            obs.measure(sw,ha,b,psi.svd(),energy);

            } //for loop over b
        
        if(obs.checkDone(sw,psi.svd(),energy)) break;
    
        } //for loop over sw
    
    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    psi.noise(orig_noise); 

    return energy;
    }


#undef Cout
#undef Endl
#undef Format

#endif
