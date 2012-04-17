//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_BASE_DMRG_WORKER_H
#define __ITENSOR_BASE_DMRG_WORKER_H

#include "DMRGOpts.h" // <-- default implementation
#include "Sweeps.h"

template <class MPSType, class DefaultOpts=DMRGOpts>
class BaseDMRGWorker
    {
    public:

    typedef typename MPSType::MPOType
    MPOType;

    BaseDMRGWorker(const Sweeps& sweeps);

    BaseDMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts);

    virtual 
    ~BaseDMRGWorker();

    const Sweeps& 
    sweeps() const { return *sweeps_; }
    void
    sweeps(const Sweeps& nswps) { sweeps_ = &nswps; }

    BaseDMRGOpts& 
    opts() const;
    void 
    opts(BaseDMRGOpts& new_opts);

    Real 
    run(const MPOType& H, MPSType& psi) 
        { return runInternal(H,psi); }

    Real 
    run(const std::vector<MPOType>& H, MPSType& psi) 
        { return runInternal(H,psi); }

    Real 
    energy() const { return getEnergy(); }

private:

    virtual Real 
    runInternal(const MPOType& H, MPSType& psi) = 0;

    virtual Real 
    runInternal(const std::vector<MPOType>& H, MPSType& psi)
        {
        Error("DMRG with a vector of MPOs not implemented for this DMRGWorker.");
        return 0;
        }

    virtual Real 
    getEnergy() const = 0;

    const Sweeps* sweeps_;

    bool own_opts_;

    BaseDMRGOpts* opts_;
    };


template <class MPSType, class DefaultOpts>
inline BaseDMRGWorker<MPSType, DefaultOpts>::
BaseDMRGWorker(const Sweeps& sweeps)
    : sweeps_(&sweeps),
      own_opts_(true),
      opts_(new DefaultOpts())
    { }

template <class MPSType, class DefaultOpts>
inline BaseDMRGWorker<MPSType, DefaultOpts>::
BaseDMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts)
    : sweeps_(&sweeps),
      own_opts_(false),
      opts_(&opts)
    { }

template <class MPSType, class DefaultOpts>
inline BaseDMRGWorker<MPSType, DefaultOpts>::
~BaseDMRGWorker()
    {
    if(own_opts_) { delete opts_; }
    }

template <class MPSType, class DefaultOpts>
inline void BaseDMRGWorker<MPSType, DefaultOpts>::
opts(BaseDMRGOpts& new_opts)
    { 
    if(own_opts_) { delete opts_; }
    opts_ = &new_opts; 
    }

template <class MPSType, class DefaultOpts>
inline BaseDMRGOpts& BaseDMRGWorker<MPSType, DefaultOpts>::
opts() const 
    { 
    return *opts_; 
    }

#endif
