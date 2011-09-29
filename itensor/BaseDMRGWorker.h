#ifndef __ITENSOR_BASE_DMRG_WORKER_H
#define __ITENSOR_BASE_DMRG_WORKER_H

#include "DMRGOpts.h" // <-- default implementation
#include "Sweeps.h"

template <class MPSType, class MPOType, 
          class DefaultDMRGOptsType=DMRGOpts>
class BaseDMRGWorker
{
public:
    BaseDMRGWorker(const Sweeps& sweeps) :
    sweeps_(sweeps),
    own_opts_(true),opts_(new DefaultDMRGOptsType())
    { }

    BaseDMRGWorker(const Sweeps& sweeps, BaseDMRGOpts& opts) :
    sweeps_(sweeps),
    own_opts_(false),opts_(opts)
    { }

    virtual ~BaseDMRGWorker()
    {
        if(own_opts_) { delete opts_; }
    }

    virtual const Sweeps& sweeps() const { return sweeps_; }

    virtual void opts(BaseDMRGOpts& new_opts)
    { 
        if(own_opts_) { delete opts_; }
        opts_ = &new_opts; 
    }
    virtual BaseDMRGOpts& opts() const { return *opts_; }

    virtual void run(const MPOType& H, MPSType& psi) = 0;
    virtual Real energy() const = 0;

private:
    const Sweeps& sweeps_;
    bool own_opts_;
    BaseDMRGOpts* opts_;
};

#endif
