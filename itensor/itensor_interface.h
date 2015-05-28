//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_INTERFACE_H
#define __ITENSOR_ITENSOR_INTERFACE_H
#include "itensor/indexset.h"
#include "itensor/detail/algs.h"

namespace itensor {

struct ITData;
template <typename DType>
struct ITDataType;

//
// IndexT - interface template for ITensor and IQTensor
//

template<typename IndexT>
class ITensorT
    {
    public:
    using index_type = IndexT;
    using indexval_type = typename index_type::indexval_type;
    using indexset_type = IndexSetT<IndexT>;
    using storage_ptr = std::shared_ptr<ITData>;
    using const_storage_ptr = std::shared_ptr<const ITData>;
    private:
    indexset_type is_;
    storage_ptr store_;
    LogNumber scale_;
    public:

    //
    // Constructors
    //

    //Default constructed tensor will evaluate to false in boolean context
    ITensorT() { }

    //Construct rank 1 tensor, all elements set to zero
    explicit
    ITensorT(const IndexT& i1) { }

    //Construct rank 2 tensor, all elements set to zero
    ITensorT(const IndexT& i1,
             const IndexT& i2) { }

    //Construct rank n tensor, all elements set to zero
    template <typename... Indices>
    ITensorT(const IndexT& i1,
             const IndexT& i2, 
             const IndexT& i3, 
             const Indices&... rest) { }

    //Construct rank 0 tensor (scalar), value set to val
    //If val.imag()==0, only Real storage will be used
    explicit
    ITensorT(Cplx val) { }

    //Construct rank n tensor, all
    //elements set to zero except the single
    //entry specified by the IndexVal/IQIndexVal args
    template <typename... IVals>
    explicit
    ITensorT(const indexval_type& iv1, 
             const IVals&... rest) { }

    operator ITensorT<Index>() const { return *this; }

    //
    // Accessor Methods
    //

    //Tensor rank (number of indices)
    int 
    r() const { return is_.r(); }

    //Access index set
    const IndexSetT<IndexT>&
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
    //of IndexVals. Will not switch storage
    //from Real to Complex unless val.imag()!=0 
    template<typename... IndexVals>
    void
    set(Cplx val, const IndexVals&... ivs) { }

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
    fill(Complex z) { return *this; }

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
    conj() { return *this; }

    ITensorT&
    dag() { return *this; }

    //Replace data with real part
    ITensorT&
    takeReal() { return *this; }

    //Replace data with imaginary part
    ITensorT&
    takeImag() { return *this; }

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
    ITensorT(IndexSetT<IndexT> iset,
             DataType&& dat,
             const LogNumber& scale = 1);

    ITensorT(IndexSetT<IndexT> iset,
             storage_ptr&& pdat,
             const LogNumber& scale = 1);

    //Provide indices from IndexSet
    explicit
    ITensorT(const IndexSetT<IndexT>& is) { }

    //Scale factor, used internally for efficient scalar ops.
    //Mostly for developer use; not necessary to explicitly involve
    //scale factors in user-level ITensor operations.
    const LogNumber&
    scale() const { return scale_; }

    LogNumber&
    scale() { return scale_; }

    const_storage_ptr
    store() const { return store_; }

    storage_ptr&
    store() { return store_; }

    const ITData&
    cstore() const { return *store_; }

    void 
    scaleTo(const LogNumber& newscale) { }

    }; // class ITensorT

template<typename IndexT>
ITensorT<IndexT>::
ITensorT(indexset_type iset,
        storage_ptr&& pdat,
        const LogNumber& scale)
    :
    is_(std::move(iset)),
    store_(std::move(pdat)),
    scale_(scale)
    { }

template<typename IndexT>
template <class DataType>
ITensorT<IndexT>::
ITensorT(indexset_type iset,
         DataType&& dat,
         const LogNumber& scale) :
    is_(std::move(iset)),
    store_(std::make_shared<ITDataType<std::decay_t<DataType>>>(std::move(dat))),
    scale_(scale)
    {
    static_assert(std::is_rvalue_reference<decltype(std::forward<DataType>(dat))>::value,
                  "Error: cannot pass lvalues to ITensorT(...,ITDataType&& dat,...) constructor");
    }

//Multiplication by real scalar
template<typename IndexT>
ITensorT<IndexT>& 
operator*=(ITensorT<IndexT>& T, Real fac)
    {
    if(fac == 0)
        {
        T.fill(0);
        return T;
        }
    T.scale() *= fac;
    return T;
    }

//Division by real scalar
template<typename IndexT>
ITensorT<IndexT>& 
operator/=(ITensorT<IndexT>& T, Real fac) { T.scale()/=fac; return T; }

//Negation
template<typename IndexT>
ITensorT<IndexT>
operator-(ITensorT<IndexT> T) { T.scale() *= -1; return T; }


template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
prime(ITensorT<IndexT> A, 
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
primeExcept(ITensorT<IndexT> A, 
            VarArgs&&... vargs)
    {
    A.primeExcept(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
noprime(ITensorT<IndexT> A, 
        VarArgs&&... vargs)
    {
    A.noprime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
mapprime(ITensorT<IndexT> A, 
         VarArgs&&... vargs)
    {
    A.mapprime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename IndexT>
bool
hasindex(const ITensorT<IndexT>& T, const typename ITensorT<IndexT>::index_type& I)
    {
    return detail::contains(T.inds(),I);
    }

template<typename IndexT>
IndexT
findtype(const ITensorT<IndexT>& T, IndexType type)
    {
    for(auto& i : T.inds())
        if(i.type()==type) return i;
    return IndexT{};
    }

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
template<typename IndexT> 
IndexT
commonIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t = All)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return IndexT{};
    }


//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
template<typename IndexT> 
IndexT
uniqueIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && !hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return IndexT{};
    }

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
          IndexType type = All)
    { 
    int tempLevel = 99999;
#ifdef DEBUG
    for(auto& I : T.inds())
        {
        if(I.primeLevel() == tempLevel) 
            {
            Print(tempLevel);
            Error("swapPrime fails if an index has primeLevel==tempLevel");
            }
        }
#endif
    T.mapprime(plev1,tempLevel,type);
    T.mapprime(plev2,plev1,type);
    T.mapprime(tempLevel,plev2,type);
    return T; 
    }

//Apply x = f(x) for each element x of T
//and return the resulting tensor
template<typename IndexT, typename F>
ITensorT<IndexT>
apply(ITensorT<IndexT> T, F&& f)
    {
    T.apply(std::forward<F>(f));
    return T;
    }

} //namespace itensor


#endif
