//
// Distributed under the ITensor Library License, Version 1.1
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDIAG_H
#define __ITENSOR_ITDIAG_H

#include "itdata.h"
#include "../simpletensor.h"

namespace itensor {

template<typename T>
class ITDiag : public ITDispatch<ITDiag<T>>
    {
    public:

    std::vector<T> data;

    template<typename... Args>
    ITDiag(Args&&... args) : data(std::forward<Args>(args)...) { }

    virtual
    ~ITDiag() { }

    };

}; //namespace itensor

#endif

