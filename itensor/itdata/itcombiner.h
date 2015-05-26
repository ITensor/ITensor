//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITCOMBINER_H
#define __ITENSOR_ITCOMBINER_H

#include <iostream>

namespace itensor {

class ITCombiner
    {
    public:

    ITCombiner() { }

    virtual
    ~ITCombiner() { }

    };

void inline
read(std::istream& s, ITCombiner& dat) { }

void inline
write(std::ostream& s, const ITCombiner& dat) { }

} //namespace itensor

#endif

