//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICERANGE_H_
#define __ITENSOR_SLICERANGE_H_

#include "itensor/tensor/range.h"

namespace itensor {

template<typename Range_, typename Perm_,
         typename RetType = Range>
RetType
permuteRange(Range_ && R,
             Perm_ const& P)
    {
    auto Rb = RangeBuilderT<RetType>(R.r());
    size_t n = 0;
    for(auto pn : P)
        {
        Rb.setIndStr(pn,R.extent(n),R.stride(n));
        ++n;
        }
    return Rb.build();
    }

template<typename RetType,
         typename Range_, typename Perm_>
RetType
permuteRangeTo(Range_ && R,
               Perm_ const& P)
    {
    return permuteRange<Range_,Perm_,RetType>(std::forward<Range_>(R),P);
    }

template<typename Range_, typename Perm_>
Range
permuteExtents(Range_ && R,
               Perm_ const& P)
    {
    auto Rb = RangeBuilder(R.r());
    size_t n = 0;
    for(auto pn : P)
        {
        Rb.setIndex(pn,R.extent(n));
        ++n;
        }
    return Rb.build();
    }

Range inline
groupIndsRange(Range const& R,
               size_t istart,
               size_t iend)
    {
    if(not isContiguous(R)) Error("groupInds requires contiguous range");
    auto ngroup = iend-istart;
    size_t nr = R.r()-ngroup+1;
    auto rb = RangeBuilder(nr);
    for(decltype(istart) j = 0; j < istart; ++j) rb.nextIndex(R.extent(j));
    auto group_ext = 1;
    for(auto j = istart; j < iend; ++j) group_ext *= R.extent(j);
    rb.nextIndex(group_ext);
    for(auto j = iend; j < size_t(R.r()); ++j) rb.nextIndex(R.extent(j));
    return rb.build();
    }

} //namespace itensor

#endif
