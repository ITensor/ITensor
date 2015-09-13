//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICERANGE_H_
#define __ITENSOR_SLICERANGE_H_

#include "itensor/matrix/range.h"

namespace itensor {

template<typename Range_, typename Perm_>
auto
permuteRange(Range_ && R,
             Perm_ const& P)
    {
    auto Rb = RangeBuilder(R.r());
    size_t n = 0;
    for(auto pn : P)
        {
        Rb.setExtStr(pn,R.extent(n),R.stride(n));
        ++n;
        }
    return Rb.build();
    }

template<typename Range_, typename Perm_>
auto
permuteExtents(Range_ && R,
               Perm_ const& P)
    {
    auto Rb = RangeBuilder(R.r());
    size_t n = 0;
    for(auto pn : P)
        {
        Rb.setExtent(pn,R.extent(n));
        ++n;
        }
    return Rb.build();
    }

} //namespace itensor

#endif
