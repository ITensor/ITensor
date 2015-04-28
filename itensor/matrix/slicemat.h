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
    return makeRef(std::forward<Mat_>(M),MRange(M.Ncols(),M.colStride(),M.Nrows(),M.rowStride()));
    }

template<typename Mat_>
auto
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
    return makeRef(std::forward<Mat_>(M),offset,subind);
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
auto
diagonal(Mat_&& M)
    {
    return makeVecRef(std::forward<Mat_>(M),std::min(M.Nrows(),M.Ncols()),M.rowStride()+M.colStride());
    }

template<typename Mat_>
auto
row(Mat_&& M, long j)
    {
    return makeVecRef(std::forward<Mat_>(M),(j-1)*M.rowStride(),M.Ncols(),M.colStride());
    }

template<typename Mat_>
auto
column(Mat_&& M, long j)
    {
    return makeVecRef(std::forward<Mat_>(M),(j-1)*M.colStride(),M.Nrows(),M.rowStride());
    }


}; //namespace itensor

#endif
