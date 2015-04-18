//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_FUNCTIONS_H
#define __ITENSOR_ITDATA_FUNCTIONS_H
#include "global.h"
#include "itdata/itdense.h"
#include "itdata/itdiag.h"
#include "itdata/itcombiner.h"
#include "itdata/iqtdata.h"
#include "indexset.h"
#include "simpletensor.h"
#include "contract.h"

namespace itensor {








//template<typename T, int size>
//struct GetPtrElt : RegisterFunc<GetPtrElt<T,size>>
//    {
//    using Inds = std::array<long,size>;
//
//    T* ptr_;
//    const Inds& inds_;
//
//    GetPtrElt(const Inds& inds)
//        : inds_(inds)
//        { }
//
//    explicit operator T*() const { return ptr_; }
//
//    template <typename V,
//              typename std::enable_if<std::is_same<V,typename std::remove_const<T>::type>::value>::type* = nullptr>
//    void
//    operator()(const ITDense<V>& d)
//        {
//        ptr_ = &(d.data.vref(d.data.ind(inds_)));
//        }
//
//    template <class D>
//    void
//    operator()(const D& d)
//        {
//        throw ITError("ITensor does not have requested element type");
//        }
//    };








}; //namespace itensor

#endif

