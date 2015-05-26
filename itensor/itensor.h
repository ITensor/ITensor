//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "itensor/itensor_interface.h"
#include "itensor/matrix/mat.h"
#include "itensor/detail/algs.h"

namespace itensor {

//
// ITensor
//
// For the ITensor class interface, see
// itensor_interface.h
// An ITensor is defined as ITensorT<Index>
//
using ITensor = ITensorT<Index>;

std::ostream& 
operator<<(std::ostream & s, const ITensor& T);

//ITensor inline
//operator*(ITensor A, const ITensor& B) { A *= B; return A; }
//ITensor inline
//operator*(const ITensor& A, ITensor&& B) { B *= A; return B; }
//ITensor inline
//operator*(ITensor T, Real fac) { T *= fac; return T; }
//ITensor inline
//operator*(Real fac, ITensor T) { T *= fac; return T; }
//ITensor inline
//operator*(ITensor T, Complex fac) { T *= fac; return T; }
//ITensor inline
//operator*(Complex fac, ITensor T) { T *= fac; return T; }
//ITensor inline
//operator/(ITensor T, Real fac) { T /= fac; return T; }
//ITensor inline
//operator/(ITensor T, Complex fac) { T /= fac; return T; }
//ITensor inline
//operator+(ITensor A, const ITensor& B) { A += B; return A; }
//ITensor inline
//operator+(const ITensor& A, ITensor&& B) { B += A; return B; }
//ITensor inline
//operator-(ITensor A, const ITensor& B) { A -= B; return A; }
//ITensor inline
//operator-(const ITensor& A, ITensor&& B) { B -= A; B *= -1; return B; }
//
//ITensor inline
//operator*(ITensor T, const IndexVal& iv) { T *= iv; return T; }
//ITensor inline
//operator*(const IndexVal& iv, const ITensor& t) { return (ITensor(iv) *= t); }

ITensor
combiner(std::vector<Index> inds, const Args& args = Global::args());

template<typename... Inds>
ITensor
combiner(const Index& i1, const Inds&... inds)
    {
    return combiner(std::vector<Index>{i1,inds...});
    }

ITensor
deltaTensor(const Index& i1, const Index& i2);

//Construct ITensor with diagonal elements set to z
//(if z is a Real or z.imag()==0 storage will be real)
template<typename... Inds>
ITensor
diagTensor(Complex z,
           const Index& i1,
           Inds&&... inds);

template<typename... Inds>
ITensor
diagTensor(Real r,
           const Index& i1,
           Inds&&... inds)
    {
    return diagTensor(Complex(r),i1,std::forward<Inds>(inds)...);
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

////
//// Define product of IndexVal iv1 = (I1,n1), iv2 = (I2,n2)
//// (I1, I2 are Index objects; n1,n2 are type int)
//// to be an ITensor T such that T(I1(n1),I2(n2)) == 1
////
//// Useful for creating MPOs
////
//ITensor inline
//operator*(const IndexVal& iv1, const IndexVal& iv2) 
//    { 
//    ITensor t(iv1); 
//    return (t *= iv2); 
//    }

//
// Define product of IndexVal iv1 = (I1,n1) with a Real "val"
// to be an ITensor T such that T(I1(n1)) == val
//
// Useful for creating MPOs
//
//ITensor inline
//operator*(const IndexVal& iv1, Real val) 
//    { 
//    ITensor res(iv1); 
//    res *= val; 
//    return res; 
//    }
//ITensor inline
//operator*(Real val, const IndexVal& iv) { return operator*(iv,val); }


template<class Tensor>
bool
hasindex(const Tensor& T, const typename Tensor::IndexT& I)
    {
    return detail::contains(T.inds(),I);
    }

ITensor
randomize(ITensor T, const Args& args = Global::args());

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
ITensor
randomTensor(const IndexSet& inds);

ITensor
matrixTensor(Mat&& M, const Index& i1, const Index& i2);

//Apply x = f(x) for each element x of T
//and return the resulting tensor
template<typename F>
ITensor
apply(ITensor T, F&& f);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//(and similar for complex case) but computed more efficiently
Real 
norm(const ITensor& T);

ITensor
conj(ITensor T);

ITensor inline
dag(const ITensor& T) { return conj(T); }

bool
isComplex(const ITensor& T);

template <class Tensor>
Tensor
realPart(Tensor T) { T.takeReal(); return T; }

template <class Tensor>
Tensor
imagPart(Tensor T) { T.takeImag(); return T; }

Real
sumels(const ITensor& T);

Complex
sumelsC(const ITensor& T);

// Read ITensor from binary input stream.
void
read(std::istream& s, ITensor& t);

// Write ITensor to binary output stream.
void 
write(std::ostream& s, const ITensor& t);


//
//Return copy of a tensor with primeLevels plev1 and plev2 swapped
//
//For example, if T has indices i,i' (like a matrix or a site
//operator) then swapPrime(T,0,1) will have indices i',i 
//i.e. the transpose of T.
//
template <class Tensor>
Tensor
swapPrime(Tensor T, int plev1, int plev2,
          IndexType type = All);

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
template<class TensorA, class TensorB> typename 
TensorA::IndexT
commonIndex(const TensorA& A, const TensorB& B, IndexType t = All);

//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
template<class TensorA, class TensorB> typename 
TensorA::IndexT
uniqueIndex(const TensorA& A, 
            const TensorB& B, 
            IndexType t);

template<class Tensor>
auto
findtype(const Tensor& T, IndexType type)
    {
    using IndexT = typename Tensor::IndexT;
    for(auto& i : T.inds())
        if(i.type()==type) return i;
    return IndexT();
    }

//
// Given Tensors which represent operator matrices
// (e.g. A(site1',site1), B(site1',site1) )
// multiply them, automatically adjusting primeLevels
// so that result is again an operator matrix C(site1',site1)
//
//              s'  t'
//  s'  t'      |   |
//  |   |       [-A-]
//  [-C-]  =    |   |
//  |   |       [-B-]
//  s   t       |   |
//              s   t
//
// (here s and t are indices of type Site)
//
template<class IndexT>
ITensorT<IndexT>
multSiteOps(ITensorT<IndexT> A, const ITensorT<IndexT>& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }



} //namespace itensor

//See file itensor.ih for template/inline method implementations
#include "itensor/itensor.ih"


#endif
