//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_INTERFACE_H
#define __ITENSOR_ITENSOR_INTERFACE_H
#include <memory>
#include "itensor/indexset.h"

namespace itensor {

struct ITData;

template<typename IndexT>
class ITensorT
    {
    public:

    using index_type = IndexT;
    using indexval_type = typename IndexT::indexval_type;
    using storage_ptr = std::shared_ptr<ITData>;

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
    ITensorT(const Index& i1, 
             const Index& i2, 
             const Index& i3, 
             const Indices&... rest);

    //Construct rank 0 tensor (scalar), value set to val
    //If val.imag()==0, only Real storage will be used
    explicit
    ITensorT(Cplx val);

    //Construct rank n tensor, all
    //elements set to zero except the single
    //entry specified by the IndexVal/IQIndexVal args
    template <typename... IVals>
    explicit
    ITensorT(const indexval_type& iv1, 
             const IVals&... rest);

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
    set(Cplx val, const IndexVals&... ivs);

    //
    // Operators
    //

    //Contracting product
    //All matching Index pairs automatically contracted
    //Cji = \sum_{k,l} Akjl * Blki
    //ITensor& 
    //operator*=(const ITensor& other);

    //// Contract with IndexVal
    //// If iv = (J,n), Index J is fixed to it's nth
    //// value and rank decreases by 1
    //// (similar to summing against a Kronecker
    //// delta tensor \delta_{J,n})
    //ITensor& 
    //operator*=(const IndexVal& iv) { return operator*=(ITensor(iv)); } 

    ////Multiplication by scalar
    //ITensor& 
    //operator*=(Real fac);
    //ITensor& 
    //operator*=(Complex z);

    ////Division by scalar
    //ITensor& 
    //operator/=(Real fac) { scale_/=fac; return *this; }
    //ITensor& 
    //operator/=(Complex z) { return operator*=(1./z); }

    ////Negation
    //ITensor
    //operator-() const { auto T = *this; T.scale_ *= -1; return T; }


    ////Tensor addition and subtraction
    ////Summands must have same Indices, in any order
    ////Cijk = Aijk + Bkij
    //ITensor& 
    //operator+=(const ITensor& other);

    //ITensor& 
    //operator-=(const ITensor& other);

    //
    // Index Prime Level Methods
    //

    template<typename... VarArgs>
    ITensorT& 
    noprime(VarArgs&&...);

    template<typename... VarArgs>
    ITensorT& 
    prime(VarArgs&&...);

    template<typename... VarArgs>
    ITensorT&
    primeExcept(VarArgs&&...);

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
    //ITensorT&
    //fill(Complex z);

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
    dag() { return conj(); }

    //Replace data with real part
    //ITensorT&
    //takeReal();

    //Replace data with imaginary part
    //ITensorT&
    //takeImag();


    private:

    //void
    //scaleOutNorm();

    //void
    //equalizeScales(ITensorT& other);

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
             const LogNumber& scale = 1) { }

    ITensorT(IndexSetT<IndexT> iset,
             storage_ptr&& pdat,
             const LogNumber& scale = 1);

    //Provide indices from IndexSet
    explicit
    ITensorT(const IndexSetT<IndexT>& is);

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
    IndexSetT<IndexT> is_;
    storage_ptr store_;
    LogNumber scale_;
    }; // class ITensorT

template<typename IndexT>
template<typename... VarArgs>
ITensorT<IndexT>& ITensorT<IndexT>::
noprime(VarArgs&&... vargs) { itensor::noprime(is_,std::forward<VarArgs>(vargs)...); return *this; }

template<typename IndexT>
template<typename... VarArgs>
ITensorT<IndexT>& ITensorT<IndexT>::
prime(VarArgs&&... vargs) { itensor::prime(is_,std::forward<VarArgs>(vargs)...); return *this; }

template<typename IndexT>
template<typename... VarArgs>
ITensorT<IndexT>& ITensorT<IndexT>::
primeExcept(VarArgs&&... vargs) { itensor::primeExcept(is_,std::forward<VarArgs>(vargs)...); return *this; }

} //namespace itensor


#endif
