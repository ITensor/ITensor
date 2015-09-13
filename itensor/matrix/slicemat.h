//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICEMAT_H_
#define __ITENSOR_SLICEMAT_H_

#include "mat.h"

namespace itensor {

template<typename Mat_>
auto
transpose(Mat_&& M)
    {
    return makeRef(std::forward<Mat_>(M),transpose(M.range()));
    }

template<typename Mat_>
auto
subMatrix(Mat_&& M,
          size_t rstart,
          size_t rstop,
          size_t cstart,
          size_t cstop)
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to subMatrix");
#ifdef DEBUG
    if(rstop > nrows(M) || rstart > rstop) throw std::runtime_error("subMatrix invalid row start and stop");
    if(cstop > ncols(M) || cstart > cstop) throw std::runtime_error("subMatrix invalid col start and stop");
#endif
    auto offset = rowStride(M)*(rstart-1)+colStride(M)*(cstart-1);
    auto subrange = MatRange(rstop-rstart+1,rowStride(M),cstop-cstart+1,colStride(M));
    return makeRef(std::forward<Mat_>(M).store()+offset,std::move(subrange));
    }

template<typename Mat_>
auto
rows(Mat_&& M,
     size_t rstart,
     size_t rstop)
    {
    return subMatrix(std::forward<Mat_>(M),rstart,rstop,1,nrows(M));
    }

template<typename Mat_>
auto
columns(Mat_&& M,
        size_t cstart,
        size_t cstop)
    {
    return subMatrix(std::forward<Mat_>(M),1,nrows(M),cstart,cstop);
    }

template<typename Mat_>
auto
diagonal(Mat_&& M)
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to diagonal(M)");
    auto drange = VecRange(std::min(nrows(M),ncols(M)),rowStride(M)+colStride(M));
    return makeRef(std::forward<Mat_>(M).store(),std::move(drange));
    }

template<typename Mat_>
auto
row(Mat_&& M, size_t j)
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to row(M,n)");
#ifdef DEBUG
    if(j < 1 || j > nrows(M)) throw std::runtime_error("invalid row index");
#endif
    auto offset = (j-1)*rowStride(M);
    return makeRef(std::forward<Mat_>(M).store()+offset,VecRange(ncols(M),colStride(M)));
    }

template<typename Mat_>
auto
column(Mat_&& M, size_t j)
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to column(M,n)");
#ifdef DEBUG
    if(j < 1 || j > ncols(M)) throw std::runtime_error("invalid column index");
#endif
    auto offset = (j-1)*colStride(M);
    return makeRef(std::forward<Mat_>(M).store()+offset,VecRange(nrows(M),rowStride(M)));
    }


} //namespace itensor

#endif
