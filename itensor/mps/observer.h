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
    measure(const Args& args = Global::args()) { }
    
    bool virtual
    checkDone(const Args& args = Global::args()) { return false; }

    virtual ~Observer() { }

    };

} //namespace itensor

#endif 
