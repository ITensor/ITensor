//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "itensor/itensor_interface.h"
#include "itensor/tensor/mat.h"
#include "itensor/tensor/slicerange.h"

namespace itensor {

//
// ITensor
//
// For the ITensor class interface, see itensor_interface.h
// For the available operators see below
//

// Contract with IndexVal
// If iv = (J,n), Index J is fixed to it's nth
// value and rank decreases by 1
// (similar to summing against a Kronecker
// delta tensor \delta_{J,n})
ITensor& 
operator*=(ITensor & T, IndexVal const& iv);
ITensor inline
operator*(ITensor T, IndexVal const& iv) { T *= iv; return T; }
ITensor
operator*(IndexVal const& iv, ITensor const& B);

ITensor
combiner(std::vector<Index> inds, Args const& args = Args::global());

template<typename... Inds>
ITensor
combiner(Index const& i1, 
         Inds const&... inds)
    {
    return combiner(std::vector<Index>{i1,inds...});
    }

Index
combinedIndex(ITensor const& C);


//Construct diagonal ITensor with diagonal 
//elements set to 1.0
template<typename... Inds>
ITensor
delta(Index const& i1,
      Inds const&... inds);

//Construct diagonal ITensor,
//diagonal elements given by container C
//(Uses elements C.begin() up to C.end())
template<typename Container, 
         typename... Inds,
         class = stdx::enable_if_t<stdx::containerOf<Real,Container>::value
                                || stdx::containerOf<Cplx,Container>::value> >
ITensor
diagTensor(Container const& C,
           Index const& i1,
           Inds&&... inds);

//
// Define product of IndexVal iv1 = (I1,n1), iv2 = (I2,n2)
// (I1, I2 are Index objects; n1,n2 are type int)
// to be an ITensor T such that T(I1(n1),I2(n2)) == 1
//
// Useful for creating MPOs
//
ITensor
operator*(IndexVal const& iv1, IndexVal const& iv2);

//
// Define product of IndexVal iv1 = (I1,n1) 
// with a scalar "fac" to be an ITensor T such that T(I1(n1)) == val
//
// Useful for creating MPOs
//
ITensor
operator*(IndexVal const& iv1, Cplx val);
ITensor inline
operator*(Cplx val, IndexVal const& iv) { return operator*(iv,val); }


template <typename... Inds>
ITensor
randomTensor(Index const& i1, Inds&&... inds)
    {
    return random(ITensor(i1,std::forward<Inds>(inds)...));
    }
template <typename... Inds>
ITensor
randomTensorC(Index const& i1, Inds&&... inds)
    {
    return random(ITensor(i1,std::forward<Inds>(inds)...),{"Complex",true});
    }

template<typename IndexT>
ITensorT<IndexT>
randomTensor(IndexSetT<IndexT> const& inds);

ITensor
matrixTensor(Matrix && M, Index const& i1, Index const& i2);
ITensor
matrixTensor(Matrix const& M, Index const& i1, Index const& i2);
ITensor
matrixTensor(CMatrix && M, Index const& i1, Index const& i2);
ITensor
matrixTensor(CMatrix const& M, Index const& i1, Index const& i2);


template<typename... Indxs>
TensorRef1
ordered(ITensor & T, Indxs&&... inds);

template<typename... Indxs>
CTensorRef1
orderedC(ITensor & T, Indxs&&... inds);

std::ostream& 
operator<<(std::ostream & s, ITensor const& T);

template<> ITensor::
ITensorT(Cplx val);

template<> inline
ITensor& ITensor::
dag() { return conj(); }


} //namespace itensor

//See file itensor.ih for template/inline method implementations
#include "itensor/itensor.ih"


#endif
