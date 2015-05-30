//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_INTERFACE_H
#define __ITENSOR_ITENSOR_INTERFACE_H
#include "itensor/detail/algs.h"
#include "itensor/itdata/dotask.h"
#include "itensor/indexset.h"

namespace itensor {

template<typename IndexT> class ITensorT;
class IQIndex;

//
// ITensorT - interface template for ITensor and IQTensor
//
// ITensor is ITensorT<Index>
// IQTensor is ITensorT<IQIndex>
//

using ITensor = ITensorT<Index>;
using IQTensor = ITensorT<IQIndex>;


struct ITData;

template<typename IndexT>
class ITensorT
    {
    public:
    using index_type = IndexT;
    using indexval_type = typename index_type::indexval_type;
    using indexset_type = IndexSetT<IndexT>;
    using storage_ptr = std::shared_ptr<ITData>;
    using const_storage_ptr = std::shared_ptr<const ITData>;
    using scale_type = LogNumber;
    private:
    indexset_type is_;
    storage_ptr store_;
    scale_type scale_;
    public:

    //
    // Constructors
    //

    //Default constructed tensor will evaluate to false in boolean context
    ITensorT() { }

    //Construct rank 1 tensor, all elements set to zero
    explicit
    ITensorT(const IndexT& i1);

    //Construct rank 2 tensor, all elements set to zero
    ITensorT(const IndexT& i1,
             const IndexT& i2);

    //Construct rank n tensor, all elements set to zero
    template <typename... Indices>
    ITensorT(const IndexT& i1,
             const IndexT& i2, 
             const IndexT& i3, 
             const Indices&... rest);

    //Construct rank 0 tensor (scalar), value set to val
    //If val.imag()==0, storage will be Real
    explicit
    ITensorT(Cplx val);

    //Construct rank n tensor, all
    //elements set to zero except the single
    //entry specified by the IndexVal/IQIndexVal args
    template <typename... IVals>
    explicit
    ITensorT(const indexval_type& iv1, 
             const IVals&... rest);

    //Automatic conversion to ITensor
    operator ITensor() const;

    //
    // Accessor Methods
    //

    //Tensor rank (number of indices)
    int 
    r() const { return is_.r(); }

    //Access index set
    const indexset_type&
    inds() const { return is_; }

    //evaluates to false if storage not allocated
    explicit operator bool() const { return bool(store_); }

    template <typename... IndexVals>
    Real
    real(IndexVals&&... ivs) const;

    template <typename... IndexVals>
    Cplx
    cplx(IndexVals&&... ivs) const;

    //Set element at location given by collection
    //of IndexVals or IQIndexVals. Will not switch storage
    //from Real to Complex unless val.imag()!=0 
    template<typename... IVals>
    void
    set(Cplx val, IVals&&... ivs);

    //
    // Index Prime Level Methods
    //

    template<typename... VarArgs>
    ITensorT& 
    noprime(VarArgs&&... vargs)
        { itensor::noprime(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensorT& 
    prime(VarArgs&&... vargs)
        { itensor::prime(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensorT&
    primeExcept(VarArgs&&... vargs)
        { itensor::primeExcept(is_,std::forward<VarArgs>(vargs)...); return *this; }

    //Change all Indices having primeLevel plevold to have primeLevel plevnew
    ITensorT& 
    mapprime(int plevold, int plevnew, IndexType type = All)
        { itensor::mapprime(is_,plevold,plevnew,type); return *this; }

    //
    // Element Transformation Methods
    //

    //Set all elements to z. If z.imag()==0
    //(such as if z is automatically converted from a Real)
    //then storage will be real only.
    ITensorT&
    fill(Complex z);

    //Call a function of the form f()->val once
    //for each element, assign result to each element.
    template <typename Func>
    ITensorT&
    generate(Func&& f);

    //Apply a function of the form f(x)->y
    //to each element x, replacing it with y
    template <typename Func>
    ITensorT&
    apply(Func&& f);

    //Apply a function of the form f(x)->void
    //to each element x.
    template <typename Func>
    const ITensorT&
    visit(Func&& f) const;

    //
    // Complex number methods
    //

    //Take complex conjugate of all elements
    ITensorT&
    conj();

    ITensorT&
    dag();

    //Replace data with real part
    ITensorT&
    takeReal();

    //Replace data with imaginary part
    ITensorT&
    takeImag();

    //
    // Operators
    //

    //Contracting product
    //All matching Index pairs automatically contracted
    //Cji = \sum_{k,l} Akjl * Blki
    ITensorT&
    operator*=(const ITensorT& other);

    //Tensor addition and subtraction
    //Summands must have same Indices, in any order
    //Cijk = Aijk + Bkij
    ITensorT& 
    operator+=(const ITensorT& other);
    ITensorT& 
    operator-=(const ITensorT& other);

    //Multiplication by real scalar
    ITensorT&
    operator*=(Real fac) { scale_ *= fac; return *this; }

    //Division by real scalar
    ITensorT&
    operator/=(Real fac) { scale_ /= fac; return *this; }

    //Multiplication by complex scalar
    ITensorT&
    operator*=(Cplx z);

    //Division by complex scalar
    ITensorT&
    operator/=(Cplx z) { return operator*=(1./z); }

    //Negation
    ITensorT&
    operator-() { scale_.negate(); return *this; }

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
    ITensorT(indexset_type iset,
             DataType&& dat,
             const scale_type& scale = 1);

    ITensorT(indexset_type iset,
             storage_ptr&& pdat,
             const scale_type& scale = 1);

    //Provide indices from IndexSet
    explicit
    ITensorT(const indexset_type& is);

    //Scale factor, used internally for efficient scalar ops.
    //Mostly for developer use; not necessary to explicitly involve
    //scale factors in user-level ITensor operations.
    const scale_type&
    scale() const { return scale_; }

    scale_type&
    scale() { return scale_; }

    const_storage_ptr
    store() const { return store_; }

    storage_ptr&
    store() { return store_; }

    void 
    scaleTo(const scale_type& newscale);

    }; // class ITensorT



template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
prime(ITensorT<IndexT> A, 
      VarArgs&&... vargs);

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
primeExcept(ITensorT<IndexT> A, 
            VarArgs&&... vargs);

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
noprime(ITensorT<IndexT> A, 
        VarArgs&&... vargs);

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
mapprime(ITensorT<IndexT> A, 
         VarArgs&&... vargs);

template<typename IndexT>
bool
hasindex(const ITensorT<IndexT>& T, const typename ITensorT<IndexT>::index_type& I);

template<typename IndexT>
IndexT
findtype(const ITensorT<IndexT>& T, IndexType type);

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
template<typename IndexT> 
IndexT
commonIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t = All);


//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
template<typename IndexT> 
IndexT
uniqueIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t);

//
//Return copy of a tensor with primeLevels plev1 and plev2 swapped
//
//For example, if T has indices i,i' (like a matrix or a site
//operator) then swapPrime(T,0,1) will have indices i',i 
//i.e. the transpose of T.
//
template <typename IndexT>
ITensorT<IndexT>
swapPrime(ITensorT<IndexT> T, 
          int plev1, 
          int plev2,
          IndexType type = All);

//Apply x = f(x) for each element x of T
//and return the resulting tensor
template<typename I, typename F>
ITensorT<I>
apply(ITensorT<I> T, F&& f);

template <class I>
ITensorT<I>
realPart(ITensorT<I> T) { T.takeReal(); return T; }

template <class I>
ITensorT<I>
imagPart(ITensorT<I> T) { T.takeImag(); return T; }

template<typename I>
bool
isComplex(const ITensorT<I>& T);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//(and similar for complex case) but computed more efficiently
template<typename I>
Real
norm(const ITensorT<I>& T);

template<typename I>
ITensorT<I>
randomize(ITensorT<I> T, const Args& args = Global::args());

template<typename I>
ITensorT<I>
conj(ITensorT<I> T);

template<typename I>
ITensorT<I>
dag(ITensorT<I> T);

template<typename I>
Real
sumels(const ITensorT<I>& t);

template<typename I>
Cplx
sumelsC(const ITensorT<I>& t);

template<typename I>
void
read(std::istream& s, ITensorT<I>& T);

template<typename I>
void
write(std::ostream& s, const ITensorT<I>& T);

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
template<class I>
ITensorT<I>
multSiteOps(ITensorT<I> A, const ITensorT<I>& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }

} //namespace itensor

#include "itensor_interface.ih"


#endif
