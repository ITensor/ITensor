//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICETEN_H_
#define __ITENSOR_SLICETEN_H_

#include "ten.h"

namespace itensor {

template<typename Ten_, typename C1, typename C2>
auto
subTensorImpl(Ten_ && T,
              C1 const& start,
              C2 const& stop)
    {
    using range_type = typename std::decay<Ten_>::type::range_type;
    using rstorage_type = typename range_type::storage_type;
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
    auto rstore = rstorage_type(r);
    auto st = start.begin();
    auto sp = stop.begin();
    for(decltype(r) j = 0; j < r; ++j, ++st, ++sp) 
        {
        offset += T.stride(j) * (*st);
        rstore[j].str = T.stride(j);
        rstore[j].ext = (*sp)-(*st);
        }
    return makeTenRef(T.data()+offset,range_type{std::move(rstore)});
    }

template<typename Ten_, typename C1, typename C2>
auto
subTensor(Ten_ && T,
          C1 const& start,
          C2 const& stop)
    {
    return subTensorImpl(std::forward<Ten_>(T),start,stop);
    }

template<typename Ten_>
auto
subTensor(Ten_ && T,
          std::initializer_list<size_t> start,
          std::initializer_list<size_t> stop)
    {
    return subTensorImpl(std::forward<Ten_>(T),start,stop);
    }


} //namespace itensor

#endif
