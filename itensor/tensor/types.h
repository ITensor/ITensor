//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TENSOR_TYPES_H
#define __ITENSOR_TENSOR_TYPES_H

#include <iostream>
#include "itensor/util/infarray.h"
#include "itensor/util/vararray.h"

namespace itensor {

using Label = InfArray<long,11ul>; //sizeof(InfArray<long,11ul>)==128
//using Label = VarArray<long,15ul>; //sizeof(VarArray<long,15ul>)==128
//using Label = VarArray<long,31ul>; //sizeof(VarArray<long,31ul>)==256

template<typename T, size_t N>
std::ostream& 
operator<<(std::ostream & s, InfArray<T,N> const& v)
    {
    if(v.empty()) return s;
    decltype(v.size()) j = 0;
    for(; 1+j < v.size(); ++j)
        {
        s << v[j] << ",";
        }
    s << v[j];
    return s;
    }

} //namespace itensor

#endif
