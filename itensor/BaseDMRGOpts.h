//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#ifndef BASE_DMRG_OPTS_H
#define BASE_DMRG_OPTS_H
#include "svdworker.h"
#include "option.h"

// virtual base class

class BaseDMRGOpts 
    {
    public:
    
    void virtual
    measure(int sw, int ha, int b, const SVDWorker& svd, Real energy,
            const Option& opt1 = Option(), const Option& opt2 = Option(),
            const Option& opt3 = Option(), const Option& opt4 = Option()) = 0;
    
    bool virtual
    checkDone(int sw, const SVDWorker& svd, Real energy, 
              const Option& opt1 = Option(), const Option& opt2 = Option()) = 0;

    virtual ~BaseDMRGOpts() { }

    };

#endif // BASE_DMRG_OPTS_H
