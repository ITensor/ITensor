//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OBSERVER_H
#define __ITENSOR_OBSERVER_H

#include "spectrum.h"

namespace itensor {

class SVDWorker;

class Observer 
    {
    public:

    Observer() { }
    
    void virtual
    measure(const OptSet& opts = Global::opts()) { }
    
    bool virtual
    checkDone(const OptSet& opts = Global::opts()) { return false; }

    virtual ~Observer() { }

    };

}; //namespace itensor

#endif 
