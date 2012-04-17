//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DMRG_WORKER_H
#define __ITENSOR_DMRG_WORKER_H
#include "BaseDMRGWorker.h"
#include "itensor.h"
#include "localmpo.h"
#include "localmposet.h"
#include "eigensolver.h"

//
// DMRGWorker
//

template <class MPSType>
class DMRGWorker : public BaseDMRGWorker<MPSType>
    {
    public:

    typedef BaseDMRGWorker<MPSType>
    Parent;

    typedef typename Parent::MPOType
    MPOType;
    
    DMRGWorker(const Sweeps& sweeps);

    DMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts);

    using Parent::sweeps;

    using Parent::opts;

    typedef typename MPSType::TensorT 
    Tensor;

    virtual 
    ~DMRGWorker() { }

    private:

    Real virtual
    runInternal(const MPOType& H, MPSType& psi);

    Real virtual 
    runInternal(const std::vector<MPOType>& H, MPSType& psi);

    virtual Real 
    getEnergy() const { return energy_; } 

    Real energy_;
   

    }; // class DMRGWorker

//
// Convenience templated wrapper functions
// so that one can call dmrg directly instead
// of using a DMRGWorker instance.
//

//DMRG with an MPO
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps)
    {
    DMRGWorker<MPSType> worker(sweeps);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with an MPO and options
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps, 
     BaseDMRGOpts& opts)
    {
    DMRGWorker<MPSType> worker(sweeps,opts);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with a set of MPOs (lazily summed)
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const std::vector<MPOType>& H, const Sweeps& sweeps)
    {
    DMRGWorker<MPSType> worker(sweeps);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with a set of MPOs and options
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const std::vector<MPOType>& H, const Sweeps& sweeps, 
     BaseDMRGOpts& opts)
    {
    DMRGWorker<MPSType> worker(sweeps,opts);
    worker.run(H,psi);
    return worker.energy();
    }


//
// DMRGWorker
//


template <class MPSType> inline
DMRGWorker<MPSType>::
DMRGWorker(const Sweeps& sweeps)
    : Parent(sweeps), energy_(0) 
    { }

template <class MPSType> inline
DMRGWorker<MPSType>::
DMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts)
    : Parent(sweeps, opts), energy_(0) 
    { }


template <class MPSType> inline
Real DMRGWorker<MPSType>::
runInternal(const MPOType& H, MPSType& psi)
    {
    typedef typename MPOType::TensorT 
    MPOTensor;

    const Real orig_cutoff = psi.cutoff(),
               orig_noise  = psi.noise();
    const int orig_minm = psi.minm(), 
              orig_maxm = psi.maxm();
    int debuglevel = (opts().quiet() ? 0 : 1);

    int N = psi.NN();
    energy_ = 0;

    psi.position(1);
    
    LocalMPO<MPOTensor> PH(H);

    Eigensolver solver;
    solver.debugLevel(debuglevel);
    
    for(int sw = 1; sw <= sweeps().nsweep(); ++sw)
        {
        psi.cutoff(sweeps().cutoff(sw)); 
        psi.minm(sweeps().minm(sw)); 
        psi.maxm(sweeps().maxm(sw));
        psi.noise(sweeps().noise(sw));
        solver.maxIter(sweeps().niter(sw));

        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            if(!opts().quiet()) 
                {
                std::cout << 
                    boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)") 
                    % sw % ha % b % (b+1) << std::endl;
                }

            PH.position(b,psi);

            Tensor phi = psi.bondTensor(b);

            energy_ = solver.davidson(PH,phi);
            
            psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH);

            if(!opts().quiet()) 
                { 
                std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d") 
                % sweeps().cutoff(sw) % sweeps().minm(sw) % sweeps().maxm(sw) << std::endl;
                std::cout << boost::format("    Trunc. err=%.1E, States kept=%s")
                % psi.svd().truncerr(b) % psi.LinkInd(b).showm() << std::endl;
                }

            opts().measure(sw,ha,b,psi.svd(),energy_);

            } //for loop over b
        
        if(opts().checkDone(sw,energy_)) break;
    
        } //for loop over sw
    
    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    psi.noise(orig_noise); 

    return energy_;
    }

template <class MPSType> inline
Real DMRGWorker<MPSType>::
runInternal(const std::vector<MPOType>& H, MPSType& psi)
    {
    typedef typename MPOType::TensorT 
    MPOTensor;

    const Real orig_cutoff = psi.cutoff(),
               orig_noise  = psi.noise();
    const int orig_minm = psi.minm(), 
              orig_maxm = psi.maxm();
    int debuglevel = (opts().quiet() ? 0 : 1);

    int N = psi.NN();
    energy_ = 0;

    psi.position(1);
    
    LocalMPOSet<MPOTensor> PH(H);

    Eigensolver solver;
    solver.debugLevel(debuglevel);
    
    for(int sw = 1; sw <= sweeps().nsweep(); ++sw)
        {
        psi.cutoff(sweeps().cutoff(sw)); 
        psi.minm(sweeps().minm(sw)); 
        psi.maxm(sweeps().maxm(sw));
        psi.noise(sweeps().noise(sw));
        solver.maxIter(sweeps().niter(sw));

        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            if(!opts().quiet()) 
                {
                std::cout << 
                    boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)") 
                    % sw % ha % b % (b+1) << std::endl;
                }

            PH.position(b,psi);

            Tensor phi = psi.bondTensor(b);

            energy_ = solver.davidson(PH,phi);
            
            psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH);

            if(!opts().quiet()) 
                { 
                std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d")
                % sweeps().cutoff(sw) % sweeps().minm(sw) % sweeps().maxm(sw) << std::endl;
                std::cout << boost::format("    Trunc. err=%.1E, States kept=%s")
                % psi.svd().truncerr(b) % psi.LinkInd(b).showm() << std::endl;
                }

            opts().measure(sw,ha,b,psi.svd(),energy_);

            } //for loop over b
        
        if(opts().checkDone(sw,energy_)) break;
    
        } //for loop over sw
    
    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    psi.noise(orig_noise); 

    return energy_;
    }


#endif // __ITENSOR_DMRG_WORKER_H
