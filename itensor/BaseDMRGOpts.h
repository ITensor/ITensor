#ifndef BASE_DMRG_OPTS_H
#define BASE_DMRG_OPTS_H
#include "mps.h"

// virtual base class

class BaseDMRGOpts {
public:
    
    virtual void measure(int sw, int ha, int b, const SVDWorker& svd, Real energy)=0;
    
    virtual bool checkDone(int sw, Real energy)=0;
};

#endif // BASE_DMRG_OPTS_H