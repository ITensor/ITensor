//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_OBSERVER_H
#define __ITENSOR_OBSERVER_H

#include "itensor/global.h"

namespace itensor {

class Observer 
    {
    public:

    Observer() { }
    
    void virtual
    measure(Args const& args = Args::global()) { }
    
    bool virtual
    checkDone(Args const& args = Args::global()) { return false; }

    virtual ~Observer() { }

    };

} //namespace itensor

#endif 
