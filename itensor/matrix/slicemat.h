//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICEMAT_H_
#define __ITENSOR_SLICEMAT_H_

#include "mat.h"

namespace itensor {

template<template<class> class Ref>
struct CheckNoRValue;

template<typename MatType>
using MSlice = std::result_of_t<CheckNoRValue<MatRefT>(MatType)>;

template<typename MatType>
using VSlice = std::result_of_t<CheckNoRValue<VecRefT>(MatType)>;

template<typename Mat_>
MSlice<Mat_>
transpose(Mat_&& M)
    {
    return makeMatRef(M.data(),MRange(M.Ncols(),M.colStride(),M.Nrows(),M.rowStride()));
    }

template<typename Mat_>
MSlice<Mat_>
subMatrix(Mat_&& M,
          long rstart,
          long rstop,
          long cstart,
          long cstop)
    {
#ifdef DEBUG
    if(rstop > M.Nrows() || rstart > rstop) throw std::runtime_error("subMatrix invalid row start and stop");
    if(cstop > M.Ncols() || cstart > cstop) throw std::runtime_error("subMatrix invalid col start and stop");
#endif
    auto offset = M.rowStride()*(rstart-1)+M.colStride()*(cstart-1);
    auto subind = MRange(rstop-rstart+1,M.rowStride(),cstop-cstart+1,M.colStride());
    return makeMatRef(M.data()+offset,subind);
    }

template<typename Mat_>
auto
rows(Mat_&& M,
     long rstart,
     long rstop)
    {
    return subMatrix(std::forward<Mat_>(M),rstart,rstop,1,M.Nrows());
    }

template<typename Mat_>
auto
columns(Mat_&& M,
        long cstart,
        long cstop)
    {
    return subMatrix(std::forward<Mat_>(M),1,M.Nrows(),cstart,cstop);
    }

template<typename Mat_>
VSlice<Mat_>
diagonal(Mat_&& M)
    {
    return makeVecRef(M.data(),std::min(M.Nrows(),M.Ncols()),M.rowStride()+M.colStride());
    }

template<typename Mat_>
VSlice<Mat_>
row(Mat_&& M, long j)
    {
    return makeVecRef(M.data()+(j-1)*M.rowStride(),M.Ncols(),M.colStride());
    }

template<typename Mat_>
VSlice<Mat_>
column(Mat_&& M, long j)
    {
    return makeVecRef(M.data()+(j-1)*M.colStride(),M.Nrows(),M.rowStride());
    }

template<typename T1, typename T2>
using transfer_const = std::conditional_t<!std::is_const<T1>::value,T2,const T2>;

template<template<class> class Ref>
struct CheckNoRValue
    {
    template<typename Mat_>
    Ref<transfer_const<Mat_,Real>>
    operator()(Mat_& m);

    template<typename D>
    Ref<D>
    operator()(Ref<D> m);
    };

}; //namespace itensor

#endif
