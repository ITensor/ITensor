//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TENSOR_TYPES_H
#define __ITENSOR_TENSOR_TYPES_H

#include <iostream>
#include "itensor/util/vararray.h"

namespace itensor {

using Label = InfArray<long,11ul>; //sizeof(InfArray<long,11ul>)==128
//using Label = VarArray<long,15ul>; //sizeof(VarArray<long,15ul>)==128
//using Label = VarArray<long,31ul>; //sizeof(VarArray<long,31ul>)==256

inline std::ostream& 
operator<<(std::ostream& s, const Label& A)
    {
    for(auto& a : A) s << a << " "; 
    s << "\n";
    return s;
    }

} //namespace itensor

#endif
