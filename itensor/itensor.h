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
inline ITensor& 
operator*=(ITensor & T, IndexVal const& iv) { return T *= ITensor(iv); } 
ITensor inline
operator*(ITensor T, IndexVal const& iv) { T *= iv; return T; }
ITensor inline
operator*(IndexVal const& iv, ITensor const& B) { ITensor A(iv); A *= B; return A; }

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
deltaTensor(const Index& i1, const Index& i2);

//Construct ITensor with diagonal elements set to z
//(if z is a Real or z.imag()==0 storage will be real)
template<typename... Inds>
ITensor
diagTensor(Cplx z,
           const Index& i1,
           Inds&&... inds);

template<typename... Inds>
ITensor
diagTensor(Real r,
           const Index& i1,
           Inds&&... inds)
    {
    return diagTensor(Cplx{r},i1,std::forward<Inds>(inds)...);
    }

//Construct diagonal ITensor,
//diagonal elements given by container C
template<typename Container, typename... Inds>
auto //->return type is ITensor
diagTensor(const Container& C,
           const Index& i1,
           Inds&&... inds) 
        //This is a "throwaway" test: we don't care about the results, just want to filter out "Container"
        //types (such as Container==int) that don't have a value_type member type
        -> typename std::conditional<std::is_same<typename Container::value_type,Real>::value,ITensor,ITensor>::type;

//
// Define product of IndexVal iv1 = (I1,n1), iv2 = (I2,n2)
// (I1, I2 are Index objects; n1,n2 are type int)
// to be an ITensor T such that T(I1(n1),I2(n2)) == 1
//
// Useful for creating MPOs
//
ITensor inline
operator*(const IndexVal& iv1, const IndexVal& iv2) 
    { 
    ITensor t(iv1); 
    return (t *= iv2); 
    }

//
// Define product of IndexVal iv1 = (I1,n1) 
// with a Real number "fac" to be an ITensor T such that T(I1(n1)) == val
//
// Useful for creating MPOs
//
ITensor inline
operator*(const IndexVal& iv1, Real val) 
    { 
    ITensor res(iv1); 
    res *= val; 
    return res; 
    }
ITensor inline
operator*(Real val, const IndexVal& iv) { return operator*(iv,val); }


template <typename... Inds>
ITensor
randomTensor(const Index& i1, Inds&&... inds)
    {
    return randomize(ITensor(i1,std::forward<Inds>(inds)...));
    }
template <typename... Inds>
ITensor
randomTensorC(const Index& i1, Inds&&... inds)
    {
    return randomize(ITensor(i1,std::forward<Inds>(inds)...),"Complex");
    }

template<typename IndexT>
ITensorT<IndexT>
randomTensor(const IndexSetT<IndexT>& inds);

ITensor
matrixTensor(Matrix&& M, const Index& i1, const Index& i2);

std::ostream& 
operator<<(std::ostream & s, const ITensor& T);



template<> ITensor::
ITensorT(const Index& i1);

template<> ITensor::
ITensorT(const Index& i1,const Index& i2);

template<> ITensor::
ITensorT(Cplx val);

template<> inline
ITensor& ITensor::
dag() { return conj(); }

} //namespace itensor

//See file itensor.ih for template/inline method implementations
#include "itensor/itensor.ih"


#endif
