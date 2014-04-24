//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OBSERVER_H
#define __ITENSOR_OBSERVER_H

#include "spectrum.h"

// virtual base class

class SVDWorker;

class Observer 
    {
    public:
    
    void virtual
    measure(int N, int sw, int ha, int b, const Spectrum& spec, Real energy,
            const OptSet& opts = Global::opts()) = 0;
    
    bool virtual
    checkDone(int sw, Real energy, 
              const OptSet& opts = Global::opts()) = 0;

    virtual ~Observer() { }

    };

#endif 
