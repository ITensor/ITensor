//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BASE_DMRG_WORKER_H
#define __ITENSOR_BASE_DMRG_WORKER_H

#include "DMRGOpts.h" // <-- default implementation
#include "Sweeps.h"

template <class MPSType, class MPOType>
class BaseDMRGWorker
    {
    public:

    BaseDMRGWorker(const Sweeps& sweeps);

    BaseDMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts);

    virtual 
    ~BaseDMRGWorker();

    const Sweeps& 
    sweeps() const { return sweeps_; }

    void 
    opts(BaseDMRGOpts& new_opts);

    BaseDMRGOpts& 
    opts() const;

    Real 
    run(const MPOType& H, MPSType& psi) 
        { return runInternal(H,psi); }

    Real 
    energy() const { return getEnergy(); }

private:

    virtual Real 
    runInternal(const MPOType& H, MPSType& psi) = 0;

    virtual Real 
    getEnergy() const = 0;

    const Sweeps& sweeps_;

    bool own_opts_;

    BaseDMRGOpts* opts_;
};


template <class MPSType, class MPOType>
inline BaseDMRGWorker<MPSType, MPOType>::
BaseDMRGWorker(const Sweeps& sweeps)
    : sweeps_(sweeps),
      own_opts_(true),
      opts_(new DMRGOpts())
    { }

template <class MPSType, class MPOType>
inline BaseDMRGWorker<MPSType, MPOType>::
BaseDMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts)
    : sweeps_(sweeps),
      own_opts_(false),
      opts_(&opts)
    { }

template <class MPSType, class MPOType>
inline BaseDMRGWorker<MPSType, MPOType>::
~BaseDMRGWorker()
    {
    if(own_opts_) { delete opts_; }
    }

template <class MPSType, class MPOType>
inline void BaseDMRGWorker<MPSType, MPOType>::
opts(BaseDMRGOpts& new_opts)
    { 
    if(own_opts_) { delete opts_; }
    opts_ = &new_opts; 
    }

template <class MPSType, class MPOType>
inline BaseDMRGOpts& BaseDMRGWorker<MPSType, MPOType>::
opts() const 
    { 
    return *opts_; 
    }

#endif
