//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICETEN_H_
#define __ITENSOR_SLICETEN_H_

#include "itensor/util/count.h"
#include "itensor/matrix/ten.h"
#include "itensor/matrix/slicerange.h"

namespace itensor {

template<typename Ten_, typename C1, typename C2>
ref_type<Ten_>
subTensor(Ten_ && T,
          C1 const& start,
          C2 const& stop);

template<typename Ten_>
ref_type<Ten_>
subIndex(Ten_ && T,
         size_t ind,
         size_t start,
         size_t stop);

//group contiguous indices
template<typename Ten_>
ref_type<Ten_>
groupInds(Ten_ && T,
          size_t istart,
          size_t iend);

//group non-contiguous indices
template<typename Ten_, typename Inds_>
Tensor
groupInds(Ten_      && T,
          Inds_ const& inds);

template<typename Ten_, typename Perm_>
ref_type<Ten_>
permute(Ten_  const& t,
        Perm_ const& P);


///
/// Implementations
/// 

template<typename Ten_, typename C1, typename C2>
ref_type<Ten_>
subTensor(Ten_ && T,
          C1 const& start,
          C2 const& stop)
    {
    static_assert(!std::is_same<Ten_&&,Tensor&&>::value,"Cannot pass temp/rvalue Tensor to subTensor");
    static_assert(!std::is_same<Ten_&&,Vector&&>::value,"Cannot pass temp/rvalue Vector to subTensor");
    static_assert(!std::is_same<Ten_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to subTensor");
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
        rb.setIndStr(j,(*sp)-(*st),T.stride(j));
        }
    return makeTenRef(makeRef(std::forward<Ten_>(T)).store()+offset,rb.build());
    }

template<typename Ten_>
ref_type<Ten_>
subIndex(Ten_ && T,
         size_t ind,
         size_t start,
         size_t stop)
    {
    static_assert(!std::is_same<Ten_&&,Tensor&&>::value,"Cannot pass temp/rvalue Tensor to subIndex");
    static_assert(!std::is_same<Ten_&&,Vector&&>::value,"Cannot pass temp/rvalue Vector to subIndex");
    static_assert(!std::is_same<Ten_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to subIndex");
#ifdef DEBUG
    if(ind >= size_t(T.r())) throw std::runtime_error("subIndex: index out of range");
#endif
    auto R = T.range();
    R[ind].ind = stop-start;
    return makeTenRef(makeRef(std::forward<Ten_>(T)).store()+T.stride(ind)*start,std::move(R));
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

template<typename Ten_>
ref_type<Ten_>
groupInds(Ten_ && T,
          size_t istart,
          size_t iend)
    {
    return makeRef(std::forward<Ten_>(T),groupIndsRange(T.range(),istart,iend));
    }

template<typename Ten_, typename Inds_>
Tensor
groupInds(Ten_      && T,
          Inds_ const& inds)
    {
    //Does permute followed by contiguous groupInds; returns a Tensor
    using value_t = decltype(inds[0]);
    auto r = T.r();
    auto P = Label(r);
    auto inds_has = [&inds](value_t j) -> long
        { 
        for(auto n : index(inds)) if(j==inds[n]) return true;
        return false;
        };
    size_t tofront = 0,
           toback = inds.size();
    for(decltype(r) j = 0; j < r; ++j)
        {
        if(inds_has(j)) P[j] = tofront++;
        else            P[j] = toback++;
        }
    auto PT = Tensor{permute(T,P)};
    if(inds.size() <= 1) return PT;
    return Tensor{std::move(PT.store()),groupIndsRange(PT.range(),0,inds.size())};
    }

template<typename Ten_, typename Perm_>
ref_type<Ten_>
permute(Ten_  && t,
        Perm_ const& P)
    {
    return makeRef(std::forward<Ten_>(t),permuteRange(t.range(),P));
    }


} //namespace itensor

#endif
