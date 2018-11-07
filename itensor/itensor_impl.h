//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_IMPL_H_
#define __ITENSOR_ITENSOR_IMPL_H_

//
// Template Method Implementations
//

namespace itensor {

//template<>
//template <typename... IVals>
//ITensor::
//ITensorT(const IndexVal& iv1, 
//         const IVals&... rest)
//  : scale_(1.)
//    {
//    const size_t size = 1+sizeof...(rest);
//    auto ivs = std::array<IndexVal,size>{{iv1,rest...}};
//    std::array<Index,size> inds;
//    for(size_t j = 0; j < size; ++j) inds[j] = ivs[j].index;
//    is_ = IndexSet(inds);
//    store_ = newITData<DenseReal>(area(is_),0.);
//    set(iv1,rest...,1.);
//    }



} //namespace itensor

#endif
