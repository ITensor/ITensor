//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "indexset.h"
#include "itdata/itdata.h"
#include "itdata/itreal.h"
#include "itdata/itcplx.h"
#include "itdata/itdiag.h"
#include "itdata/itcombiner.h"
#include "itdata/iqtdata.h"
#include "matrix/mat.h"

namespace itensor {

//
// ITensor
//
class ITensor
    {
    public:

    using IndexT = Index;
    using IndexValT = IndexVal;
    using storage_ptr = PData;

    //
    // Constructors
    //

    //Construct Null ITensor, ITensor will evaluate to false in boolean context
    ITensor() { }

    //Construct rank 1 ITensor, all elements set to zero
    explicit
    ITensor(const Index& i1);

    //Construct rank 2 ITensor, all elements set to zero
    ITensor(const Index& i1,
            const Index& i2);

    //Construct rank n ITensor, all elements set to zero
    template <typename... Indices>
    ITensor(const Index& i1, 
            const Index& i2, 
            const Index& i3, 
            const Indices&... rest);

    //Construct rank 0 ITensor (scalar), value set to val
    //If val.imag()==0, only Real storage will be used
    explicit
    ITensor(Complex val);

    //Construct rank n ITensor, all
    //elements set to zero except the single
    //entry specified by the IndexVal args
    template <typename... IVals>
    explicit
    ITensor(const IndexVal& iv1, 
            const IVals&... rest);

    //
    // Accessor Methods
    //

    //Rank of this ITensor (number of indices)
    int 
    r() const { return is_.r(); }

    //Access index set
    const IndexSet&
    inds() const { return is_; }

    //evaluate to false if ITensor is default constructed
    explicit operator bool() const { return bool(store_); }

    template <typename... IndexVals>
    Real
    real(IndexVals&&... ivs) const;

    template <typename... IndexVals>
    Complex
    cplx(IndexVals&&... ivs) const;

    //Set element at location given by collection
    //of IndexVals. Will not switch storage
    //from Real to Complex unless val.imag()!=0 
    template<typename... IndexVals>
    void
    set(Complex val, const IndexVals&... ivs);

    //
    // Operators
    //

    //Contracting product
    //All matching Index pairs automatically contracted
    //Cji = \sum_{k,l} Akjl * Blki
    ITensor& 
    operator*=(const ITensor& other);

    // Contract with IndexVal
    // If iv = (J,n), Index J is fixed to it's nth
    // value and rank decreases by 1
    // (similar to summing against a Kronecker
    // delta tensor \delta_{J,n})
    ITensor& 
    operator*=(const IndexVal& iv) { return operator*=(ITensor(iv)); } 

    //Multiplication by scalar
    ITensor& 
    operator*=(Real fac);
    ITensor& 
    operator*=(Complex z);

    //Division by scalar
    ITensor& 
    operator/=(Real fac) { scale_/=fac; return *this; }
    ITensor& 
    operator/=(Complex z) { return operator*=(1./z); }

    //Negation
    ITensor
    operator-() const { auto T = *this; T.scale_ *= -1; return T; }


    //Tensor addition and subtraction
    //Summands must have same Indices, in any order
    //Cijk = Aijk + Bkij
    ITensor& 
    operator+=(const ITensor& other);

    ITensor& 
    operator-=(const ITensor& other);

    //
    // Index Prime Level Methods
    //

    //Set primeLevel of Indices to zero
    ITensor& 
    noprime(IndexType type = All) { itensor::noprime(is_,type); return *this; }

    //Set primeLevel of Index I to zero
    template<typename... Inds>
    ITensor& 
    noprime(const Index& I1, Inds&&... inds);

    //Increase primeLevel of Indices by 1 (or optional amount inc)
    ITensor& 
    prime(int inc = 1) { prime(All,inc); return *this; }

    ITensor& 
    prime(IndexType type, int inc = 1) { itensor::prime(is_,type,inc); return *this; }

    template<typename... Inds>
    ITensor& 
    prime(const Index& I1, Inds&&... inds);

    template<typename... IVals>
    ITensor& 
    prime(const IndexVal& iv, IVals&&... ivals);

    template<typename... Inds>
    ITensor&
    primeExcept(const Index& I1, Inds&&... inds);

    template<typename... ITs>
    ITensor&
    primeExcept(IndexType it, ITs&&... etc);

    //Change all Indices having primeLevel plevold to have primeLevel plevnew
    ITensor& 
    mapprime(int plevold, int plevnew, IndexType type = All)
        { itensor::mapprime(is_,plevold,plevnew,type); return *this; }

    //
    // Element Transformation Methods
    //

    //Set all elements to z. If z.imag()==0
    //(such as if z is automatically converted from a Real)
    //then storage will be real only.
    ITensor&
    fill(Complex z);

    //Call a function of the form f()->val once
    //for each element, assign result to each element.
    template <typename Func>
    ITensor&
    generate(Func&& f);

    //Apply a function of the form f(x)->y
    //to each element x, replacing it with y
    template <typename Func>
    ITensor&
    apply(Func&& f);

    //Apply a function of the form f(x)->void
    //to each element x.
    template <typename Func>
    const ITensor&
    visit(Func&& f) const;

    //
    // Complex number methods
    //

    //Take complex conjugate of all elements
    ITensor&
    conj();

    ITensor&
    dag() { return conj(); }

    //Replace data with real part
    ITensor&
    takeReal();

    //Replace data with imaginary part
    ITensor&
    takeImag();


    private:

    void
    scaleOutNorm();

    void
    equalizeScales(ITensor& other);

    public:

    //
    // Developer / advanced methods
    //
    // The following methods should not
    // be needed for most user code.
    //

    //Construct by explicitly providing data object
    //DataType should be a subclass of ITData
    template <class DataType>
    ITensor(IndexSet iset,
            DataType&& dat,
            const LogNumber& scale = 1);

    ITensor(IndexSet iset,
            storage_ptr&& pdat,
            const LogNumber& scale = 1);

    //Provide indices from IndexSet
    explicit
    ITensor(const IndexSet& is);

    //Scale factor, used internally for efficient scalar ops.
    //Mostly for developer use; not necessary to explicitly involve
    //scale factors in user-level ITensor operations.
    const LogNumber&
    scale() const { return scale_; }

    const ITData&
    data() const { return *store_; }

    storage_ptr&
    pdata() { return store_; }

    void 
    scaleTo(const LogNumber& newscale);

    private:
    IndexSet is_;
    storage_ptr store_;
    LogNumber scale_;
    }; // class ITensor

std::ostream& 
operator<<(std::ostream & s, const ITensor& T);

// Read ITensor from binary input stream.
void
read(std::istream& s, ITensor& t);

// Write ITensor to binary output stream.
void 
write(std::ostream& s, const ITensor& t);

ITensor inline
operator*(ITensor A, const ITensor& B) { A *= B; return A; }
ITensor inline
operator*(const ITensor& A, ITensor&& B) { B *= A; return B; }
ITensor inline
operator*(ITensor T, Real fac) { T *= fac; return T; }
ITensor inline
operator*(Real fac, ITensor T) { T *= fac; return T; }
ITensor inline
operator*(ITensor T, Complex fac) { T *= fac; return T; }
ITensor inline
operator*(Complex fac, ITensor T) { T *= fac; return T; }
ITensor inline
operator/(ITensor T, Real fac) { T /= fac; return T; }
ITensor inline
operator/(ITensor T, Complex fac) { T /= fac; return T; }
ITensor inline
operator+(ITensor A, const ITensor& B) { A += B; return A; }
ITensor inline
operator+(const ITensor& A, ITensor&& B) { B += A; return B; }
ITensor inline
operator-(ITensor A, const ITensor& B) { A -= B; return A; }
ITensor inline
operator-(const ITensor& A, ITensor&& B) { B -= A; B *= -1; return B; }

ITensor inline
operator*(ITensor T, const IndexVal& iv) { T *= iv; return T; }
ITensor inline
operator*(const IndexVal& iv, const ITensor& t) { return (ITensor(iv) *= t); }

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
// Define product of IndexVal iv1 = (I1,n1) with a Real "val"
// to be an ITensor T such that T(I1(n1)) == val
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

//Return copy of ITensor A with primeLevel of indices
//I1,I2,... increased (by optional amount inc if passed
//as last argument, e.g. prime(A,i1,i2,i3,i4,5); )
template <class Tensor, typename... Inds>
Tensor
prime(Tensor A, 
      const typename Tensor::IndexT& I1, 
      Inds&&... inds);

//Return copy of ITensor A with primeLevel of indices
//iv1.index,iv2.index,... changed by iv1.val, iv2.val,...
template <class Tensor, typename... IVals>
Tensor
prime(Tensor A, 
      const typename Tensor::IndexValT& iv,
      IVals&&... ivals);

template<class Tensor, typename... Inds>
Tensor 
primeExcept(Tensor A,
            const typename Tensor::IndexT& I1, 
            Inds&&... inds);

template<class Tensor, typename... ITs>
Tensor 
primeExcept(Tensor A,
            IndexType it,
            ITs&&... etc);

//Return copy of ITensor with primeLevel of Index I set to zero
template <class Tensor, class IndexT, typename... Inds>
Tensor
noprime(Tensor A,
        const IndexT& I, 
        Inds&&... inds);

//template <class Tensor>
//Tensor
//noprime(Tensor A);

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
ITensor inline
randomTensor(const IndexSet& inds)
    {
    return randomize(ITensor(inds));
    }

ITensor
matrixTensor(Mat&& M, const Index& i1, const Index& i2);

//template <typename... Indices>
//ITensor
//tieIndex(const ITensor& T,
//         const Index& t0,
//         const Index& t1,
//         const Indices&... rest);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//[and similar for complex case] but computed much more efficiently.
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
template<class Tensor>
Tensor
multSiteOps(Tensor A, const Tensor& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }

} //namespace itensor

//See file itensor.ih for template/inline method implementations
#include "itensor.ih"


#endif
