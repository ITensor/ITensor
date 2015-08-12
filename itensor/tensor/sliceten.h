//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICETEN_H_
#define __ITENSOR_SLICETEN_H_

#include "itensor/tensor/ten.h"
#include "itensor/tensor/slicerange.h"

namespace itensor {

template<typename Ten_, typename C1, typename C2>
auto
subTensor(Ten_ && T,
          C1 const& start,
          C2 const& stop);

template<typename T, typename R, typename Perm_>
auto
permute(TenRef<T,R> const& t,
        Perm_ const& P);

template<typename T, typename R, typename Perm_>
auto
permute(Ten<T,R> const& t,
        Perm_ const& P);


///
/// Implementations
/// 

template<typename Ten_, typename C1, typename C2>
auto
subTensor(Ten_ && T,
          C1 const& start,
          C2 const& stop)
    {
    using range_type = decltype(T.range());
    using stop_type = decltype(*stop.begin());
    auto r = T.r();
#ifdef DEBUG
    if(r != decltype(r)(start.size())) throw std::runtime_error("subTensor: wrong size of start");
    if(r != decltype(r)(stop.size()))  throw std::runtime_error("subTensor: wrong size of stop");
    auto st_ = start.begin();
    auto sp_ = stop.begin();
    for(decltype(r) j = 0; j < r; ++j, ++st_, ++sp_)
        {
        if(*sp_ > stop_type(T.extent(j))) throw std::runtime_error("subTensor: stop value too large");
        if(*st_ >= *sp_)    throw std::runtime_error("subTensor: start value >= stop value");
        }
#endif
    size_t offset = 0;
    auto rb = RangeBuilderT<range_type>(r);
    auto st = start.begin();
    auto sp = stop.begin();
    for(decltype(r) j = 0; j < r; ++j, ++st, ++sp) 
        {
        offset += T.stride(j) * (*st);
        rb.setExtStr(j,(*sp)-(*st),T.stride(j));
        }
    return makeTenRef(T.data()+offset,rb.build());
    }

template<typename Ten_>
auto
subIndex(Ten_ && T,
         size_t ind,
         size_t start,
         size_t stop)
    {
#ifdef DEBUG
    if(ind >= size_t(T.r())) throw std::runtime_error("subIndex: index out of range");
#endif
    auto R = T.range();
    R[ind].ext = stop-start;
    return makeTenRef(T.data()+T.stride(ind)*start,std::move(R));
    }

template<typename Ten_>
auto
groupInds(Ten_ && T,
          size_t istart,
          size_t iend)
    {
    if(not isContiguous(T.range())) Error("groupInds requires contiguous range");
    auto ngroup = iend-istart;
    size_t nr = T.r()-ngroup+1;
    auto rb = RangeBuilder(nr);
    for(size_t j = 0; j < istart; ++j) rb.nextExtent(T.extent(j));
    auto group_ext = 1;
    for(size_t j = istart; j < iend; ++j) group_ext *= T.extent(j);
    rb.nextExtent(group_ext);
    for(size_t j = iend; j < size_t(T.r()); ++j) rb.nextExtent(T.extent(j));
    return makeTenRef(T.data(),rb.build());
    }

template<typename T, typename R, typename Perm_>
auto
permute(TenRef<T,R> const& t,
        Perm_ const& P)
    {
    return makeTenRef(t.data(),permuteRange(t.range(),P));
    }

template<typename T, typename R, typename Perm_>
auto
permute(Ten<T,R> const& t,
        Perm_ const& P)
    {
    return permute(makeRef(t),P);
    }

} //namespace itensor

#endif
