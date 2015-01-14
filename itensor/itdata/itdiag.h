//
// Distributed under the ITensor Library License, Version 1.2
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
    using size_type = typename std::vector<T>::size_type;

    T val = 0;
    std::vector<T> data;

    template<typename InputIterator>
    ITDiag(InputIterator b, InputIterator e)
        :
        data(b,e)
        { }

    ITDiag(size_t size, T val)
        :
        data(size,val)
        { }

    ITDiag(T t) 
        : val(t) 
        { }

    virtual
    ~ITDiag() { }

    bool
    allSame() const { return data.empty(); }

    };

}; //namespace itensor

#endif

