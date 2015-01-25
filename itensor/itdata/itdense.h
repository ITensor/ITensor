//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDENSE_H
#define __ITENSOR_ITDENSE_H

#include "itdata.h"

namespace itensor {

//
// TODO: replace std::vector storage with
//       storage type only holding data ptr
//       and size, maybe use in simpletensor too
//

template<typename T>
class ITDense : public ITDispatch<ITDense<T>>
    {
    public:

    std::vector<T> data;

    ITDense() { }

    ITDense(size_t size) : data(size) { }

    ITDense(size_t size, T val) : data(size,val) { }

    template<typename InputIterator>
    ITDense(InputIterator b, InputIterator e) : data(b,e) { }

    virtual
    ~ITDense() { }

    };

}; //namespace itensor

#endif

