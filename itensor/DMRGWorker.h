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
#include "localmpo_mps.h"
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
    
    DMRGWorker(const Sweeps& sweeps,
               const OptSet& opts = Global::opts());

    DMRGWorker(const Sweeps& sweeps, Observer& obs,
               const OptSet& opts = Global::opts());

    using Parent::sweeps;

    using Parent::observer;

    typedef typename MPSType::TensorT 
    Tensor;

    virtual 
    ~DMRGWorker() { }

    Real
    weight() const { return weight_; }
    void
    weight(Real val) { weight_ = val; }

    private:

    /////////////
    //
    // Data Members

    Real energy_;
    bool quiet_;
    Real weight_;

    //
    /////////////

    void
    parseOptions(const OptSet& opts);

    Real virtual
    runInternal(const MPOType& H, MPSType& psi);

    Real virtual 
    runInternal(const std::vector<MPOType>& H, MPSType& psi);

    Real virtual 
    runInternal(const MPOType& H, const std::vector<MPSType> psis, MPSType& psi);

    Real virtual 
    getEnergy() const { return energy_; } 

    }; // class DMRGWorker

//
// Convenience templated wrapper functions
// so that one can call dmrg directly instead
// of using a DMRGWorker instance.
//

//DMRG with an MPO
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps,
     const OptSet& opts = Global::opts())
    {
    DMRGWorker<MPSType> worker(sweeps,opts);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with an MPO and a custom Observer
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const MPOType& H, const Sweeps& sweeps, 
     Observer& obs,
     const OptSet& opts = Global::opts())
    {
    DMRGWorker<MPSType> worker(sweeps,obs,opts);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with a set of MPOs (lazily summed)
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const std::vector<MPOType>& H, const Sweeps& sweeps,
     const OptSet& opts = Global::opts())
    {
    DMRGWorker<MPSType> worker(sweeps,opts);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with a set of MPOs and a custom Observer
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, const std::vector<MPOType>& H, const Sweeps& sweeps, 
     Observer& obs,
     const OptSet& opts = Global::opts())
    {
    DMRGWorker<MPSType> worker(sweeps,obs,opts);
    worker.run(H,psi);
    return worker.energy();
    }

//DMRG with a single Hamiltonian MPO and a set of 
//MPS to orthogonalize against
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, 
     const MPOType& H, const std::vector<MPSType>& psis, 
     const Sweeps& sweeps, 
     const OptSet& opts = Global::opts())
    {
    DMRGWorker<MPSType> worker(sweeps,opts);
    worker.run(H,psis,psi);
    return worker.energy();
    }

//DMRG with a single Hamiltonian MPO and a set of 
//MPS to orthogonalize against, as well as a custom Observer
template <class MPSType, class MPOType>
Real inline
dmrg(MPSType& psi, 
     const MPOType& H, const std::vector<MPSType>& psis, 
     const Sweeps& sweeps, Observer& obs, 
     const OptSet& opts = Global::opts())
    {
    DMRGWorker<MPSType> worker(sweeps,obs,opts);
    worker.run(H,psis,psi);
    return worker.energy();
    }


//
// DMRGWorker
//


template <class MPSType> inline
DMRGWorker<MPSType>::
DMRGWorker(const Sweeps& sweeps,
           const OptSet& opts)
    : 
    Parent(sweeps), 
    energy_(0),
    quiet_(false),
    weight_(1)
    { 
    parseOptions(opts);
    }

template <class MPSType> inline
DMRGWorker<MPSType>::
DMRGWorker(const Sweeps& sweeps, Observer& obs,
           const OptSet& opts)
    : 
    Parent(sweeps, obs), 
    energy_(0),
    quiet_(false),
    weight_(1)
    { 
    parseOptions(opts);
    }

template <class MPSType> inline
void DMRGWorker<MPSType>::
parseOptions(const OptSet& opts)
    {
    quiet_ = opts.getBool("Quiet",false);
    weight_ = opts.getReal("Weight",1);
    }


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
    int debuglevel = (quiet_ ? 0 : 1);

    int N = psi.N();
    energy_ = NAN;

    psi.position(1);
    
    LocalMPO<MPOTensor> PH(H);

    Eigensolver solver;
    solver.debugLevel(debuglevel);

    const Opt doNorm = DoNormalize(true);
    
    for(int sw = 1; sw <= sweeps().nsweep(); ++sw)
        {
        psi.cutoff(sweeps().cutoff(sw)); 
        psi.minm(sweeps().minm(sw)); 
        psi.maxm(sweeps().maxm(sw));
        psi.noise(sweeps().noise(sw));
        solver.maxIter(sweeps().niter(sw));

        if(!PH.doWrite() 
           && Global::opts().defined("WriteM")
           && sweeps().maxm(sw) >= Global::opts().getInt("WriteM"))
            {
            std::string write_dir = Global::opts().getString("WriteDir","./");

            if(!quiet_)
                std::cout << "\nTurning on write to disk, write_dir = " << write_dir << std::endl;

            psi.doWrite(true);
            PH.doWrite(true);
            }

        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            if(!quiet_)
                {
                std::cout << 
                    boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)") 
                    % sw % ha % b % (b+1) << std::endl;
                }

            PH.position(b,psi);

            Tensor phi = psi.bondTensor(b);

            energy_ = solver.davidson(PH,phi);
            
            psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,doNorm);

            if(!quiet_)
                { 
                std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d") 
                % sweeps().cutoff(sw) % sweeps().minm(sw) % sweeps().maxm(sw) << std::endl;
                std::cout << boost::format("    Trunc. err=%.1E, States kept=%s")
                % psi.svd().truncerr(b) % psi.LinkInd(b).showm() << std::endl;
                }

            observer().measure(sw,ha,b,psi.svd(),energy_);

            } //for loop over b
        
        if(observer().checkDone(sw,psi.svd(),energy_)) break;
    
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
    int debuglevel = (quiet_ ? 0 : 1);

    int N = psi.N();
    energy_ = 0;

    psi.position(1);
    
    LocalMPOSet<MPOTensor> PH(H);

    Eigensolver solver;
    solver.debugLevel(debuglevel);

    const Opt doNorm = DoNormalize(true);
    
    for(int sw = 1; sw <= sweeps().nsweep(); ++sw)
        {
        psi.cutoff(sweeps().cutoff(sw)); 
        psi.minm(sweeps().minm(sw)); 
        psi.maxm(sweeps().maxm(sw));
        psi.noise(sweeps().noise(sw));
        solver.maxIter(sweeps().niter(sw));

        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            if(!quiet_)
                {
                std::cout << 
                    boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)") 
                    % sw % ha % b % (b+1) << std::endl;
                }

            PH.position(b,psi);

            Tensor phi = psi.bondTensor(b);

            energy_ = solver.davidson(PH,phi);
            
            psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,doNorm);

            if(!quiet_)
                { 
                std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d")
                % sweeps().cutoff(sw) % sweeps().minm(sw) % sweeps().maxm(sw) << std::endl;
                std::cout << boost::format("    Trunc. err=%.1E, States kept=%s")
                % psi.svd().truncerr(b) % psi.LinkInd(b).showm() << std::endl;
                }

            observer().measure(sw,ha,b,psi.svd(),energy_);

            } //for loop over b
        
        if(observer().checkDone(sw,psi.svd(),energy_)) break;
    
        } //for loop over sw
    
    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    psi.noise(orig_noise); 

    return energy_;
    }

template <class MPSType> inline
Real DMRGWorker<MPSType>::
runInternal(const MPOType& H, const std::vector<MPSType> psis, MPSType& psi)
    {
    typedef typename MPOType::TensorT 
    MPOTensor;

    const Real orig_cutoff = psi.cutoff(),
               orig_noise  = psi.noise();
    const int orig_minm = psi.minm(), 
              orig_maxm = psi.maxm();
    int debuglevel = (quiet_ ? 0 : 1);

    int N = psi.N();
    energy_ = 0;

    psi.position(1);
    
    LocalMPO_MPS<MPOTensor> PH(H,psis);
    PH.weight(this->weight_);

    Eigensolver solver;
    solver.debugLevel(debuglevel);

    const Opt doNorm = DoNormalize(true);
    
    for(int sw = 1; sw <= sweeps().nsweep(); ++sw)
        {
        psi.cutoff(sweeps().cutoff(sw)); 
        psi.minm(sweeps().minm(sw)); 
        psi.maxm(sweeps().maxm(sw));
        psi.noise(sweeps().noise(sw));
        solver.maxIter(sweeps().niter(sw));

        if(!PH.doWrite() 
           && Global::opts().defined("WriteM")
           && sweeps().maxm(sw) >= Global::opts().getInt("WriteM"))
            {
            std::string write_dir = Global::opts().getString("WriteDir","./");

            if(!quiet_)
                std::cout << "\nTurning on write to disk, write_dir = " << write_dir << std::endl;

            psi.doWrite(true);
            PH.doWrite(true);
            }

        for(int b = 1, ha = 1; ha != 3; sweepnext(b,ha,N))
            {
            if(!quiet_)
                {
                std::cout << 
                    boost::format("Sweep=%d, HS=%d, Bond=(%d,%d)") 
                    % sw % ha % b % (b+1) << std::endl;
                }

            PH.position(b,psi);

            Tensor phi = psi.bondTensor(b);

            energy_ = solver.davidson(PH,phi);
            
            psi.svdBond(b,phi,(ha==1?Fromleft:Fromright),PH,doNorm);

            if(!quiet_)
                { 
                std::cout << boost::format("    Truncated to Cutoff=%.1E, Min_m=%d, Max_m=%d")
                % sweeps().cutoff(sw) % sweeps().minm(sw) % sweeps().maxm(sw) << std::endl;
                std::cout << boost::format("    Trunc. err=%.1E, States kept=%s")
                % psi.svd().truncerr(b) % psi.LinkInd(b).showm() << std::endl;
                }

            observer().measure(sw,ha,b,psi.svd(),energy_);

            } //for loop over b
        
        if(observer().checkDone(sw,psi.svd(),energy_)) break;
    
        } //for loop over sw
    
    psi.cutoff(orig_cutoff); 
    psi.minm(orig_minm); 
    psi.maxm(orig_maxm);
    psi.noise(orig_noise); 

    return energy_;
    }

#endif // __ITENSOR_DMRG_WORKER_H
