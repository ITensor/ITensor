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

class ITensor;

//
// ITensor
//

struct ITData;

class ITensor
    {
    public:
    using index_type = Index;
    using indexval_type = IndexVal;
    using indexset_type = IndexSet;
    using range_type = RangeT<Index>;
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
    ITensor() { }

    //Construct rank n tensor, all elements set to zero
    //Usage: ITensor(i1,i2,i3,...)
    template <typename... index_types>
    ITensor(index_type  const& i1,
             index_types const&... i2etc);

    explicit
    ITensor(std::vector<index_type> const& inds);

    template<size_t N> 
    explicit
    ITensor(std::array<index_type,N> const& inds);

    ITensor(std::initializer_list<index_type> inds);

    //Construct rank 0 tensor (scalar), value set to val
    //If val.imag()==0, storage will be Real
    explicit
    ITensor(Cplx val);

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

    template <typename IV, typename... IVs>
    auto
    cplx(IV const& iv1, IVs&&... ivs) const
         -> stdx::if_compiles_return<Cplx,decltype(iv1.index),decltype(iv1.val)>;

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
    ITensor& 
    noprime(VarArgs&&... vargs)
        { itensor::noprime(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensor& 
    prime(VarArgs&&... vargs)
        { itensor::prime(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensor& 
    primeLevel(VarArgs&&... vargs)
        { itensor::primeLevel(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensor&
    primeExcept(VarArgs&&... vargs)
        { itensor::primeExcept(is_,std::forward<VarArgs>(vargs)...); return *this; }

    //Change all Indices having primeLevel plevold to have primeLevel plevnew
    ITensor& 
    mapprime(int plevold, int plevnew, IndexType type = All)
        { itensor::mapprime(is_,plevold,plevnew,type); return *this; }

    template<typename... VarArgs>
    ITensor& 
    mapprime(VarArgs&&... vargs)
        { itensor::mapprime(is_,std::forward<VarArgs>(vargs)...); return *this; }

    template<typename... VarArgs>
    ITensor& 
    sim(VarArgs&&... vargs)
        { itensor::sim(is_,std::forward<VarArgs>(vargs)...); return *this; }

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
    dag();

    //Replace data with real part
    ITensor&
    takeReal();

    //Replace data with imaginary part
    ITensor&
    takeImag();

    //
    // Operators
    //

    //Contracting product
    //All matching Index pairs automatically contracted
    //Cji = \sum_{k,l} Akjl * Blki
    ITensor&
    operator*=(ITensor const& other);

    //Tensor addition and subtraction
    //Summands must have same Indices, in any order
    //Cijk = Aijk + Bkij
    ITensor& 
    operator+=(ITensor const& other);
    ITensor& 
    operator-=(ITensor const& other);

#ifdef USESCALE
    //Multiplication by real scalar
    ITensor&
    operator*=(Real fac) { scale_ *= fac; return *this; }

    //Division by real scalar
    ITensor&
    operator/=(Real fac) { scale_ /= fac; return *this; }
#else
    //Multiplication by real scalar
    ITensor&
    operator*=(Real fac);

    //Division by real scalar
    ITensor&
    operator/=(Real fac);
#endif

    //Multiplication by complex scalar
    ITensor&
    operator*=(Cplx z);

    //Division by complex scalar
    ITensor&
    operator/=(Cplx z) { return operator*=(1./z); }

    //Negation
    ITensor
    operator-() const;

    //Non-contracting product
    //All matching Index pairs automatically merged
    //Ciik = Aij * Bjk
    ITensor&
    operator/=(ITensor const& other);

    //template<typename... Indxs>
    //ITensor&
    //order(index_type const& ind1, Indxs const&... inds);

    template<typename... Indxs>
    auto 
    order(index_type const& ind1, Indxs const&... inds)
    -> stdx::enable_if_t<not stdx::and_<std::is_same<index_type, Indxs>...>::value,ITensor&>;

    template <typename... Indxs>
    auto 
    order(index_type const& ind1, Indxs const&... inds)
        -> stdx::enable_if_t<stdx::and_<std::is_same<index_type, Indxs>...>::value,ITensor&>;

    template<typename... Indxs>
    ITensor&
    order(std::string const& dots, Indxs const&... inds);

    ITensor&
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
    ITensor(indexset_type iset,
             StorageType&& store,
             scale_type const& scale = LogNum{1.});

    ITensor(indexset_type iset,
             storage_ptr&& pstore,
             scale_type const& scale = LogNum{1.});

    //Provide indices from IndexSet
    explicit
    ITensor(indexset_type const& is);


    storage_ptr&
    store() { return store_; }

    const_storage_ptr
    store() const { return const_storage_ptr(store_); }


    void
    swap(ITensor & other);
    
    
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

    }; // class ITensor

//
// ITensor special constructor functions
//

// Makes a tensor with element specified by IndexVals
// set to 1.0, all other elements zero
template <typename IVal, typename... IVals>
ITensor
setElt(IVal  const& iv1, 
       IVals const&... rest);


//
// ITensor prime level functions
//
template<typename... VarArgs>
ITensor
prime(ITensor A, 
      VarArgs&&... vargs);

template<typename... VarArgs>
ITensor
primeLevel(ITensor A, 
           VarArgs&&... vargs);

template<typename... VarArgs>
ITensor
primeExcept(ITensor A, 
            VarArgs&&... vargs);

template<typename... VarArgs>
ITensor
noprime(ITensor A, 
        VarArgs&&... vargs);

template<typename... VarArgs>
ITensor
mapprime(ITensor A, 
         VarArgs&&... vargs);

template<typename... VarArgs>
ITensor
sim(ITensor A, 
    VarArgs&&... vargs);

bool
hasindex(ITensor const& T, Index const& I);

Index
findtype(ITensor const& T, IndexType type);

template<typename Cond>
Index
findindex(ITensor const& T, Cond && cond);

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
Index
commonIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t = All);


//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
Index
uniqueIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t);

//
//Return copy of a tensor with primeLevels plev1 and plev2 swapped
//
//For example, if T has indices i,i' (like a matrix or a site
//operator) then swapPrime(T,0,1) will have indices i',i 
//i.e. the transpose of T.
//
ITensor
swapPrime(ITensor T, 
          int plev1, 
          int plev2,
          IndexType type = All);

//Apply x = f(x) for each element x of T
//and return the resulting tensor
template<typename F>
ITensor
apply(ITensor T, F&& f);

ITensor inline
realPart(ITensor T) { T.takeReal(); return T; }

ITensor inline
imagPart(ITensor T) { T.takeImag(); return T; }

bool
isComplex(ITensor const& T);

bool
isReal(ITensor const& T);

//return number of indices of T
//(same as order)
long
rank(ITensor const& T);

//return number of indices of T
//(same as rank)
long
ord(ITensor const& T);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt((T*T).real()) 
//(and similar for complex case) but computed more efficiently
Real
norm(ITensor const& T);

void
randomize(ITensor & T, Args const& args = Args::global());

ITensor
random(ITensor T, Args const& args = Args::global());

ITensor
conj(ITensor T);

ITensor
dag(ITensor T);

Real
sumels(ITensor const& t);

Cplx
sumelsC(ITensor const& t);

template<typename... Inds>
ITensor
reindex(ITensor const& cT, 
        Index o1, Index n1, 
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
ITensor inline
multSiteOps(ITensor A, ITensor const& B);

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

ITensor
randomTensor(IndexSet const& inds);

ITensor
matrixTensor(Matrix && M, Index const& i1, Index const& i2);
ITensor
matrixTensor(Matrix const& M, Index const& i1, Index const& i2);
ITensor
matrixTensor(CMatrix && M, Index const& i1, Index const& i2);
ITensor
matrixTensor(CMatrix const& M, Index const& i1, Index const& i2);


//template<typename... Indxs>
//TensorRef1
//ordered(ITensor & T, Indxs&&... inds);

//template<typename... Indxs>
//CTensorRef1
//orderedC(ITensor & T, Indxs&&... inds);

std::ostream& 
operator<<(std::ostream & s, ITensor const& T);

} //namespace itensor

#include "itensor_interface_impl.h"


#endif
