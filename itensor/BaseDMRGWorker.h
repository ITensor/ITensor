//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BASE_DMRG_WORKER_H
#define __ITENSOR_BASE_DMRG_WORKER_H

#include "DMRGObserver.h" // <-- default implementation
#include "Sweeps.h"

template <class MPSType, class DefaultObserver=DMRGObserver>
class BaseDMRGWorker
    {
    public:

    typedef typename MPSType::MPOType
    MPOType;

    BaseDMRGWorker(const Sweeps& sweeps);

    BaseDMRGWorker(const Sweeps& sweeps, Observer& obs);

    virtual 
    ~BaseDMRGWorker();

    const Sweeps& 
    sweeps() const { return *sweeps_; }
    void
    sweeps(const Sweeps& nswps) { sweeps_ = &nswps; }

    Observer& 
    observer() const;
    void 
    observer(Observer& new_obs);

    Real 
    run(const MPOType& H, MPSType& psi) 
        { return runInternal(H,psi); }

    Real 
    run(const std::vector<MPOType>& H, MPSType& psi) 
        { return runInternal(H,psi); }

    Real 
    run(const MPOType& H, const std::vector<MPSType>& psis, MPSType& psi) 
        { return runInternal(H,psis,psi); }

    Real 
    energy() const { return getEnergy(); }

private:

    virtual Real 
    runInternal(const MPOType& H, MPSType& psi)
        {
        Error("DMRG for single wavefunction and Hamiltonian not implemented for this DMRGWorker.");
        return 0;
        }

    virtual Real 
    runInternal(const std::vector<MPOType>& H, MPSType& psi)
        {
        Error("DMRG with a vector of MPOs not implemented for this DMRGWorker.");
        return 0;
        }

    virtual Real 
    runInternal(const MPOType& H, const std::vector<MPSType> psis, MPSType& psi)
        {
        Error("DMRG orthogonalized against a vector of MPS not implemented for this DMRGWorker.");
        return 0;
        }

    virtual Real 
    getEnergy() const = 0;

    const Sweeps* sweeps_;

    bool own_obs_;

    Observer* obs_;
    };


template <class MPSType, class DefaultObserver>
inline BaseDMRGWorker<MPSType, DefaultObserver>::
BaseDMRGWorker(const Sweeps& sweeps)
    : sweeps_(&sweeps),
      own_obs_(true),
      obs_(new DefaultObserver())
    { }

template <class MPSType, class DefaultObserver>
inline BaseDMRGWorker<MPSType, DefaultObserver>::
BaseDMRGWorker(const Sweeps& sweeps, Observer& obs)
    : sweeps_(&sweeps),
      own_obs_(false),
      obs_(&obs)
    { }

template <class MPSType, class DefaultObserver>
inline BaseDMRGWorker<MPSType, DefaultObserver>::
~BaseDMRGWorker()
    {
    if(own_obs_) { delete obs_; }
    }

template <class MPSType, class DefaultObserver>
inline void BaseDMRGWorker<MPSType, DefaultObserver>::
observer(Observer& new_obs)
    { 
    if(own_obs_) { delete obs_; }
    obs_ = &new_obs; 
    }

template <class MPSType, class DefaultObserver>
inline Observer& BaseDMRGWorker<MPSType, DefaultObserver>::
observer() const 
    { 
    return *obs_; 
    }

#endif
