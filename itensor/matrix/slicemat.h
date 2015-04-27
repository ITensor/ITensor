//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICEMAT_H_
#define __ITENSOR_SLICEMAT_H_

#include "mat.h"

namespace itensor {

template<typename MatT>
auto
transpose(MatT& M)
    {
    auto& i = M.ind();
    return makeMatRef(M,MRange(i.cn,i.cs,i.rn,i.rs));
    }


}; //namespace itensor

#endif
