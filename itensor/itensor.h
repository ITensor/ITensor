//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "itensor/itensor_interface.h"
#include "itensor/tensor/mat.h"

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
combiner(std::vector<Index> inds, Args const& args = Global::args());

template<typename... Inds>
ITensor
combiner(Index const& i1, 
         Inds const&... inds)
    {
    return combiner(std::vector<Index>{i1,inds...});
    }

ITensor
delta(Index const& i1, Index const& i2);

//Construct ITensor with diagonal elements set to z
//(if z is a Real or z.imag()==0 storage will be real)
template<typename... Inds>
ITensor
diag(Cplx z,
     Index const& i1,
     Inds&&... inds);

template<typename... Inds>
ITensor
diag(Real r,
     Index const& i1,
     Inds&&... inds)
    {
    return diag(Cplx{r},i1,std::forward<Inds>(inds)...);
    }

//Construct diagonal ITensor,
//diagonal elements given by container C
//(Uses elements C.begin() up to C.end())
template<typename Container, 
         typename... Inds,
         class = stdx::enable_if_t<stdx::containerOf<Real,Container>::value
                                || stdx::containerOf<Cplx,Container>::value> >
ITensor
diag(Container const& C,
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
// with a Real number "fac" to be an ITensor T such that T(I1(n1)) == val
//
// Useful for creating MPOs
//
ITensor
operator*(IndexVal const& iv1, Real val);
ITensor inline
operator*(Real val, IndexVal const& iv) { return operator*(iv,val); }


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
matrixTensor(Matrix&& M, Index const& i1, Index const& i2);

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
