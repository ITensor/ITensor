//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_INTERFACE_H
#define __ITENSOR_ITENSOR_INTERFACE_H
#include "itensor/detail/algs.h"
#include "itensor/itdata/applyfunc.h"
#include "itensor/indexset.h"
#define REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/storage_types.h"

namespace itensor {

template<typename index_type> 
class ITensorT;

class IQIndex;

//
// ITensorT - interface template for ITensor and IQTensor
//
// ITensor  is ITensorT<Index>
// IQTensor is ITensorT<IQIndex>
//

using ITensor  = ITensorT<Index>;
using IQTensor = ITensorT<IQIndex>;

struct ITData;

template<typename index_type_>
class ITensorT
    {
    public:
    using index_type = index_type_;
    using indexval_type = typename index_type::indexval_type;
    using indexset_type = IndexSetT<index_type>;
    using range_type = RangeT<index_type>;
    using size_type = typename range_type::size_type;
    using storage_ptr = PData;
    using const_storage_ptr = CPData;
    using scale_type = LogNum;
    private:
    indexset_type is_;
    mutable storage_ptr store_;
    IF_USESCALE(scale_type scale_;)
    public:

    //
    // Constructors
    //

    //Default constructed tensor will evaluate to false in boolean context
    ITensorT() { }

    //Construct rank n tensor, all elements set to zero
    //Usage: ITensor(i1,i2,i3,...)
    template <typename... index_types>
    ITensorT(index_type  const& i1,
             index_types const&... i2etc);

    explicit
    ITensorT(std::vector<index_type> const& inds);

    template<size_t N> 
    explicit
    ITensorT(std::array<index_type,N> const& inds);

    ITensorT(std::initializer_list<index_type> inds);

    //Construct rank 0 tensor (scalar), value set to val
    //If val.imag()==0, storage will be Real
    explicit
    ITensorT(Cplx val);

    //Automatic conversion IQTensor -> ITensor
    operator ITensor() const;

    //
    // Accessor Methods
    //

    //Tensor rank (number of indices)
    int 
    r() const { return is_.r(); }

    //Access index set
    indexset_type const&
    inds() const { return is_; }

    //Access index
    index_type const&
    index(size_type I) const { return is_.index(I); }

    //evaluates to false if default constructed
    explicit operator bool() const { return bool(is_) || bool(store_); }

    template <typename... IndexVals>
    Real
    real(IndexVals&&... ivs) const;

    Cplx
    cplx(std::vector<indexval_type> const& ivs) const;

    template <typename IV, typename... IVs>
    auto
    cplx(IV const& iv1, IVs&&... ivs) const
         -> stdx::if_compiles_return<Cplx,decltype(iv1.index),decltype(iv1.val)>;

    template <typename Int>
    auto
    cplx(std::vector<Int> const& ints) const
        -> stdx::enable_if_t<std::is_integral<Int>::value,Cplx>;

    template <typename Int, typename... Ints>
    auto
    cplx(Int iv1, Ints... ivs) const
        -> stdx::enable_if_t<std::is_integral<Int>::value && stdx::and_<std::is_integral<Ints>...>::value,Cplx>;

    Cplx
    cplx() const;

    //Set element at location given by collection
    //of IndexVals or IQIndexVals. Will not switch storage
    //from Real to Complex unless val.imag()!=0 
    template<typename IV, typename... VArgs>
    auto
    set(IV const& iv1, VArgs&&... ivs)
        -> stdx::if_compiles_return<void,decltype(iv1.index),decltype(iv1.val)>;

    template<typename Int, typename... VArgs>
    auto
    set(Int iv1, VArgs&&... ivs)
        -> stdx::enable_if_t<std::is_integral<Int>::value,void>;

    void
    set(Cplx val);

    void
    set(std::vector<indexval_type> const& ivs, Cplx val);

    void
    set(std::vector<int> const& ivs, Cplx val);

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
    primeLevel(VarArgs&&... vargs)
        { itensor::primeLevel(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensorT&
    primeExcept(VarArgs&&... vargs)
        { itensor::primeExcept(is_,std::forward<VarArgs>(vargs)...); return *this; }

    //Change all Indices having primeLevel plevold to have primeLevel plevnew
    ITensorT& 
    mapprime(int plevold, int plevnew, IndexType type = All)
        { itensor::mapprime(is_,plevold,plevnew,type); return *this; }

    template<typename... VarArgs>
    ITensorT& 
    mapprime(VarArgs&&... vargs)
        { itensor::mapprime(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensorT& 
    sim(VarArgs&&... vargs)
        { itensor::sim(is_,std::forward<VarArgs>(vargs)...); return *this; }

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
    operator*=(ITensorT const& other);

    //Tensor addition and subtraction
    //Summands must have same Indices, in any order
    //Cijk = Aijk + Bkij
    ITensorT& 
    operator+=(ITensorT const& other);
    ITensorT& 
    operator-=(ITensorT const& other);

#ifdef USESCALE
    //Multiplication by real scalar
    ITensorT&
    operator*=(Real fac) { scale_ *= fac; return *this; }

    //Division by real scalar
    ITensorT&
    operator/=(Real fac) { scale_ /= fac; return *this; }
#else
    //Multiplication by real scalar
    ITensorT&
    operator*=(Real fac);

    //Division by real scalar
    ITensorT&
    operator/=(Real fac);
#endif

    //Multiplication by complex scalar
    ITensorT&
    operator*=(Cplx z);

    //Division by complex scalar
    ITensorT&
    operator/=(Cplx z) { return operator*=(1./z); }

    //Negation
    ITensorT
    operator-() const;

    //Non-contracting product
    //All matching Index pairs automatically merged
    //Ciik = Aij * Bjk
    ITensorT&
    operator/=(ITensorT const& other);

    //template<typename... Indxs>
    //ITensorT&
    //order(index_type const& ind1, Indxs const&... inds);

    template<typename... Indxs>
    auto 
    order(index_type const& ind1, Indxs const&... inds)
    -> stdx::enable_if_t<not stdx::and_<std::is_same<index_type, Indxs>...>::value,ITensorT&>;

    template <typename... Indxs>
    auto 
    order(index_type const& ind1, Indxs const&... inds)
        -> stdx::enable_if_t<stdx::and_<std::is_same<index_type, Indxs>...>::value,ITensorT&>;

    template<typename... Indxs>
    ITensorT&
    order(std::string const& dots, Indxs const&... inds);

    ITensorT&
    order(indexset_type const& iset);

    //
    // Read from and write to streams
    //

    void
    read(std::istream& s);

    void
    write(std::ostream& s) const;


    //
    // Developer / advanced methods
    //
    // The following methods should not
    // be needed for most user code.
    //

    template <class StorageType>
    ITensorT(indexset_type iset,
             StorageType&& store,
             scale_type const& scale = LogNum{1.});

    ITensorT(indexset_type iset,
             storage_ptr&& pstore,
             scale_type const& scale = LogNum{1.});

    //Provide indices from IndexSet
    explicit
    ITensorT(indexset_type const& is);


    storage_ptr&
    store() { return store_; }

    const_storage_ptr
    store() const { return const_storage_ptr(store_); }


    void
    swap(ITensorT & other);
    
    
#ifdef USESCALE

    scale_type const&
    scale() const { return scale_; }

    scale_type&
    scale() { return scale_; }
    
    void 
    scaleTo(scale_type const& newscale);
    
    void 
    scaleTo(Real newscale);

#else //not using scale, default case:

    scale_type
    scale() const { return scale_type(1.); }

    //scale_type&
    //scale() { return scale_; }
    
    void 
    scaleTo(scale_type const& newscale) { }
    
    void 
    scaleTo(Real newscale) { }

#endif

    }; // class ITensorT

//
// ITensorT special constructor functions
//

// Makes a tensor with element specified by IndexVals/IQIndexVals
// set to 1.0, all other elements zero
template <typename IVal, typename... IVals>
//ITensorT<typename IVal::index_type>
ITensorT<typename std::common_type<IVal,IVals...>::type::index_type>
setElt(IVal  const& iv1, 
       IVals const&... rest);


//
// ITensorT prime level functions
//
template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
prime(ITensorT<IndexT> A, 
      VarArgs&&... vargs);

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
primeLevel(ITensorT<IndexT> A, 
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

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
sim(ITensorT<IndexT> A, 
    VarArgs&&... vargs);

template<typename IndexT>
bool
hasindex(const ITensorT<IndexT>& T, const typename ITensorT<IndexT>::index_type& I);

template<typename IndexT>
IndexT
findtype(const ITensorT<IndexT>& T, IndexType type);

template<typename IndexT,
         typename Cond>
IndexT
findindex(ITensorT<IndexT> const& T, Cond && cond);

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
uniqueIndex(ITensorT<IndexT> const& A, 
            ITensorT<IndexT> const& B, 
            IndexType t);

//Find index of tensor A (of optional type t) 
//which is NOT shared by the tensors "Ts"
template<typename IndexT, typename... Tensors> 
IndexT
uniqueIndex(ITensorT<IndexT> const& A, 
            ITensorT<IndexT> const& T1,
            ITensorT<IndexT> const& T2,
            Tensors const&... Ts);

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
isComplex(ITensorT<I> const& T);

template<typename I>
bool
isReal(ITensorT<I> const& T);

//return number of indices of T
//(same as order)
template<typename I>
long
rank(ITensorT<I> const& T);

//return number of indices of T
//(same as rank)
template<typename I>
long
ord(ITensorT<I> const& T);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//(and similar for complex case) but computed more efficiently
template<typename I>
Real
norm(ITensorT<I> const& T);

template<typename I>
void
randomize(ITensorT<I> & T, Args const& args = Args::global());

template<typename I>
ITensorT<I>
random(ITensorT<I> T, Args const& args = Args::global());

template<typename I>
ITensorT<I>
conj(ITensorT<I> T);

template<typename I>
ITensorT<I>
dag(ITensorT<I> T);

template<typename I>
Real
sumels(ITensorT<I> const& t);

template<typename I>
Cplx
sumelsC(ITensorT<I> const& t);

template<typename IndexT, typename... Inds>
ITensorT<IndexT>
reindex(ITensorT<IndexT> const& cT, 
        IndexT o1, IndexT n1, 
        Inds... inds);


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
multSiteOps(ITensorT<I> A, ITensorT<I> const& B);

template<class I>
detail::IndexValIter<I>
iterInds(ITensorT<I> const& T)
    {
    return iterInds(T.inds());
    }

std::ostream& 
operator<<(std::ostream & s, ITensor const& T);

std::ostream& 
operator<<(std::ostream & s, IQTensor const& T);

} //namespace itensor

#include "itensor_interface_impl.h"


#endif
