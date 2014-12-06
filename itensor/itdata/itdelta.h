//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDELTA_H
#define __ITENSOR_ITDELTA_H

#include "itdata.h"

namespace itensor {

//
// ITDelta storage represents
// sparse storage where all diagonal
// elements are the same (== z)
// and the rest are zero
//
// Useful for replacing indices,
// grouping indices together,
// and "ungrouping" indices 
// (replacing one by many)
//
// 

class ITDelta : public ITDispatch<ITDelta>
    {
    public:

    Complex z;

    ITDelta() : z(1.) { }
    ITDelta(Real r) : z(r) { }
    ITDelta(Complex z_) : z(z_) { }

    virtual
    ~ITDelta() { }

    };

}; //namespace itensor

#endif

