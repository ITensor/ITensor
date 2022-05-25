//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "itensor/detail/algs.h"
#include "itensor/itdata/applyfunc.h"
#include "itensor/indexset.h"
#define REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/storage_types.h"
#include "itensor/tensor/mat.h"

namespace itensor {


class ITensor
    {
    public:
    using range_type = RangeT<Index>;
    using size_type = typename range_type::size_type;
    using storage_ptr = PData;
    using const_storage_ptr = CPData;
    using scale_type = LogNum;
    private:
    IndexSet is_;
    mutable storage_ptr store_;
    IF_USESCALE(scale_type scale_;)
    public:

    //
    // Constructors
    //

    //Default constructed tensor will evaluate to false in boolean context
    ITensor() { }

    explicit
    ITensor(IndexSet const& inds);

    explicit
    ITensor(std::initializer_list<Index> inds);

    //Construct n-index tensor, all elements set to zero
    //Usage: ITensor(i1,i2,i3,...)
    template <typename... Inds>
    explicit
    ITensor(Index const& i1,
            Inds const&... inds);

    //template<size_t N> 
    //explicit
    //ITensor(std::array<Index,N> const& inds);

    ITensor(QN q, IndexSet const& inds);

    template <typename... Inds>
    explicit
    ITensor(QN q, Index const& i1,
            Inds const&... inds);

    //Construct order 0 tensor (scalar), value set to val
    //If val.imag()==0, storage will be Real
    explicit
    ITensor(Cplx val);

    //
    // Accessor Methods
    //

    //Tensor order (number of indices)
    int 
    order() const { return is_.order(); }

    //Access index set
    IndexSet const&
    inds() const { return is_; }

    //Access index
    Index const&
    index(size_type I) const { return is_(I); }

    //evaluates to false if default constructed
    explicit operator bool() const { return bool(is_) || bool(store_); }

    template <typename... IndexVals>
    Real
    elt(IndexVals&&... ivs) const;

    template <typename IV, typename... IVs>
    auto
    eltC(IV const& iv1, IVs&&... ivs) const
         -> stdx::if_compiles_return<Cplx,decltype(iv1.index),decltype(iv1.val)>;

    template <typename Int, typename... Ints>
    auto
    eltC(Int iv1, Ints... ivs) const
        -> stdx::enable_if_t<std::is_integral<Int>::value 
                          && stdx::and_<std::is_integral<Ints>...>::value,Cplx>;

    Cplx
    eltC() const;

    Cplx
    eltC(std::vector<IndexVal> const& ivs) const;

    template<typename Int>
    auto
    eltC(std::vector<Int> const& ints) const
        -> stdx::enable_if_t<std::is_integral<Int>::value,Cplx>;

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
    set(std::vector<IndexVal> const& ivs, Cplx val);

    void
    set(std::vector<int> const& ivs, Cplx val);

    ITensor&
    randomize(Args const& args = Args::global());

    ITensor&
    replaceInds(IndexSet const& is1,
                IndexSet const& is2);

    ITensor&
    swapInds(IndexSet const& is1,
             IndexSet const& is2);

    //
    // Index Tag Methods
    //

    ITensor&
    setTags(TagSet const& ts,
            IndexSet const& is)
        { is_.setTags(ts,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    setTags(VarArgs&&... vargs)
        { is_.setTags(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    noTags(IndexSet const& is)
        { is_.noTags(is); return *this; }

    template<typename... VarArgs>
    ITensor&
    noTags(VarArgs&&... vargs)
        { is_.noTags(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    addTags(TagSet const& ts,
            IndexSet const& is)
        { is_.addTags(ts,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    addTags(VarArgs&&... vargs)
        { is_.addTags(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    removeTags(TagSet const& ts,
               IndexSet const& is)
        { is_.removeTags(ts,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    removeTags(VarArgs&&... vargs)
        { is_.removeTags(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    replaceTags(TagSet const& ts1,
                TagSet const& ts2,
                IndexSet const& is)
        { is_.replaceTags(ts1,ts2,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    replaceTags(VarArgs&&... vargs)
        { is_.replaceTags(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    swapTags(TagSet const& ts1,
             TagSet const& ts2,
             IndexSet const& is)
        { is_.swapTags(ts1,ts2,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    swapTags(VarArgs&&... vargs)
        { is_.swapTags(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    setPrime(int plev,
             IndexSet const& is)
        { is_.setPrime(plev,is); return *this; }

    template<typename... VarArgs>
    ITensor& 
    setPrime(VarArgs&&... vargs)
        { is_.setPrime(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    mapPrime(int plevold,
             int plevnew,
             IndexSet const& is)
        { is_.mapPrime(plevold,plevnew,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    mapPrime(VarArgs&&... vargs)
        { is_.mapPrime(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    swapPrime(int plevold,
              int plevnew,
              IndexSet const& is)
        { is_.swapPrime(plevold,plevnew,is); return *this; }

    template<typename... VarArgs>
    ITensor&
    swapPrime(VarArgs&&... vargs)
        { is_.swapPrime(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    noPrime(IndexSet const& is)
        { is_.noPrime(is); return *this; }

    template<typename... VarArgs>
    ITensor& 
    noPrime(VarArgs&&... vargs)
        { is_.noPrime(std::forward<VarArgs>(vargs)...); return *this; }

    ITensor&
    prime(int plev,
          IndexSet const& is)
        { is_.prime(plev,is); return *this; }

    ITensor&
    prime(IndexSet const& is)
        { is_.prime(is); return *this; }

    template<typename... VarArgs>
    ITensor& 
    prime(VarArgs&&... vargs)
        { is_.prime(std::forward<VarArgs>(vargs)...); return *this; }

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

    ITensor&
    fixBlockDeficient();

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

    //Force storage to be complex-valued
    ITensor&
    makeCplx();

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

    ITensor&
    permute(IndexSet const& iset);

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
    explicit
    ITensor(IndexSet iset,
            StorageType&& store,
            scale_type const& scale = LogNum{1.});

    explicit
    ITensor(IndexSet iset,
            storage_ptr&& pstore,
            scale_type const& scale = LogNum{1.});

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

    //
    // Deprecated methods
    //

    int 
    r() const { return this->order(); }

    template <typename... IndexVals>
    Real
    real(IndexVals&&... ivs) const;

    template <typename... IndexVals>
    Cplx
    cplx(IndexVals&&... ivs) const;

    template<typename... VarArgs>
    ITensor& 
    noprime(VarArgs&&... vargs)
        { Error(".noprime() is deprecated, use .noPrime() instead"); return *this; }

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

// Get ITensor values
template <typename... VarArgs>
Real
elt(ITensor A,
    VarArgs&&... vargs);

template <typename... VarArgs>
Cplx
eltC(ITensor A,
     VarArgs&&... vargs);

//
//  Templated version that we can call from template code. Simply forwards to elt or eltc based on T.
//  Real and Complex specializations are in itensor_impl.h
//
template <typename T>
T
eltT(ITensor A);     



// Get IndexSet
IndexSet const& 
inds(ITensor const& A);

long
minDim(ITensor const& A);

long
maxDim(ITensor const& A);

// Get a vector of IndexSets from a vector of ITensors
std::vector<IndexSet> 
inds(std::vector<ITensor> const& A);

size_t
nnzblocks(ITensor const& A);

long
nnz(ITensor const& A);

// Get Index at a certain position
// in the ITensor's IndexSet, using 1-based indexing
Index const& 
index(ITensor const& A, RangeT<Index>::size_type I);

ITensor
operator*(ITensor A, ITensor const& B);
ITensor
operator*(ITensor const& A, ITensor&& B);
ITensor
operator*(ITensor T, Real fac);
ITensor
operator*(Real fac, ITensor T);
ITensor
operator*(ITensor T, Complex fac);
ITensor
operator*(Complex fac, ITensor T);
ITensor
operator/(ITensor T, Real fac);
ITensor
operator/(ITensor T, Complex fac);
ITensor
operator+(ITensor A, ITensor const& B);
ITensor
operator+(ITensor const& A, ITensor&& B);
ITensor
operator-(ITensor A, ITensor const& B);
ITensor
operator-(ITensor const& A, ITensor&& B);
ITensor
operator/(ITensor A, ITensor const& B);
ITensor
operator/(ITensor const& A, ITensor && B);

// Partial direct sum of ITensors A and B
// over the specified indices
std::tuple<ITensor,IndexSet>
directSum(ITensor const& A, ITensor const& B,
          IndexSet const& I, IndexSet const& J,
          Args const& args = Args::global());

std::tuple<ITensor,Index>
directSum(ITensor const& A, ITensor const& B,
          Index const& i, Index const& j,
          Args const& args = Args::global());

//
// ITensor tag functions
//

// Define versions that explicitly take
// IndexSet for matching, so that other containers
// of Index that can be converted to IndexSet
// can be used (i.e. std::vector<Index>,
// std::initializer_list<Index>, etc.)
ITensor
setTags(ITensor A,
        TagSet const& ts,
        IndexSet const& is);

template<typename... VarArgs>
ITensor
setTags(ITensor A,
        VarArgs&&... vargs)
    {
    A.setTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
noTags(ITensor A,
       IndexSet const& is);

template<typename... VarArgs>
ITensor
noTags(ITensor A,
       VarArgs&&... vargs)
    {
    A.noTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
addTags(ITensor A,
        TagSet const& ts,
        IndexSet const& is);

template<typename... VarArgs>
ITensor
addTags(ITensor A,
        VarArgs&&... vargs)
    {
    A.addTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
removeTags(ITensor A,
           TagSet const& ts,
           IndexSet const& is);

template<typename... VarArgs>
ITensor
removeTags(ITensor A,
           VarArgs&&... vargs)
    {
    A.removeTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
replaceTags(ITensor A,
            TagSet const& ts1,
            TagSet const& ts2,
            IndexSet const& is);

template<typename... VarArgs>
ITensor
replaceTags(ITensor A,
            VarArgs&&... vargs)
    {
    A.replaceTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
swapTags(ITensor A,
         TagSet const& ts1,
         TagSet const& ts2,
         IndexSet const& is);

template<typename... VarArgs>
ITensor
swapTags(ITensor A,
         VarArgs&&... vargs)
    {
    A.swapTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
prime(ITensor A,
      int plev,
      IndexSet const& is);

ITensor
prime(ITensor A,
      IndexSet const& is);

template<typename... VarArgs>
ITensor
prime(ITensor A,
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
setPrime(ITensor A,
         int plev,
         IndexSet const& is);

template<typename... VarArgs>
ITensor
setPrime(ITensor A,
         VarArgs&&... vargs)
    {
    A.setPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
mapPrime(ITensor A,
         int plevold,
         int plevnew,
         IndexSet const& is);

template<typename... VarArgs>
ITensor
mapPrime(ITensor A,
         VarArgs&&... vargs)
    {
    A.mapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
swapPrime(ITensor A,
          int plevold,
          int plevnew,
          IndexSet const& is);

template<typename... VarArgs>
ITensor
swapPrime(ITensor A,
          VarArgs&&... vargs)
    {
    A.swapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

ITensor
noPrime(ITensor A,
        IndexSet const& is);

template<typename... VarArgs>
ITensor
noPrime(ITensor A,
        VarArgs&&... vargs)
    {
    A.noPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

bool
hasIndex(ITensor const& T,
         Index const& imatch);

bool
hasInds(ITensor const& T,
        IndexSet const& ismatch);

Arrow
dir(ITensor const& T,
    Index const& i);

template<typename... Inds>
bool
hasInds(ITensor const& T,
        Index const& i1, Inds&&... inds)
  {
  return hasInds(T,IndexSet(i1,std::forward<Inds>(inds)...));
  }

Index
findIndex(ITensor const& T,
          TagSet const& tsmatch);

IndexSet
findInds(ITensor const& T,
         TagSet const& tsmatch);

IndexSet
commonInds(ITensor const& A,
           ITensor const& B);
IndexSet
commonInds(ITensor const& A,
           ITensor const& B,
           TagSet const& tsmatch);

//Find index of tensor A (optionally having tags ts)
//which is shared with tensor B
Index
commonIndex(ITensor const& A, 
            ITensor const& B);
Index
commonIndex(ITensor const& A, 
            ITensor const& B, 
            TagSet const& tsmatch);

//Find index of tensor A
//which is NOT shared by tensor(s) B
IndexSet
uniqueInds(ITensor const& A,
           ITensor const& B);
IndexSet
uniqueInds(ITensor const& A,
           std::vector<ITensor> const& B);
IndexSet
uniqueInds(ITensor const& A,
           std::initializer_list<ITensor> B);

Index
uniqueIndex(ITensor const& A, 
            ITensor const& B);
Index
uniqueIndex(ITensor const& A, 
            ITensor const& B,
            TagSet const& tsmatch);
Index
uniqueIndex(ITensor const& A, 
            std::vector<ITensor> const& B);
Index
uniqueIndex(ITensor const& A, 
            std::vector<ITensor> const& B,
            TagSet const& tsmatch);
Index
uniqueIndex(ITensor const& A,
            std::initializer_list<ITensor> B);
Index
uniqueIndex(ITensor const& A,
            std::initializer_list<ITensor> B,
            TagSet const& tsmatch);

//permute function which returns a new ITensor 
ITensor
permute(ITensor A,
        IndexSet const& is);

template<typename... Inds>
ITensor
permute(ITensor A,
        Index const& i1,
        Inds&&... inds)
    {
    return permute(A, IndexSet(i1, std::forward<Inds>(inds)...));
    }

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
long
order(ITensor const& T);

//Compute the norm of an ITensor.
//Thinking of elements as a vector, equivalent to sqrt(v*v).
//Result is equivalent to sqrt(elt(T*T)) 
//(and similar for complex case) but computed more efficiently
Real
norm(ITensor const& T);

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

ITensor
replaceInds(ITensor T,
            IndexSet const& is1,
            IndexSet const& is2);

ITensor
swapInds(ITensor T,
         IndexSet const& is1,
         IndexSet const& is2);

detail::IndexValIter
iterInds(ITensor const& T);

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
ITensor
multSiteOps(ITensor A, ITensor const& B);

//
// Special ITensor constructors
//

std::tuple<ITensor,Index>
combiner(IndexSet const& inds, Args const& args = Args::global());

template<typename... Inds>
std::tuple<ITensor,Index>
combiner(Index const& i1, 
         Inds&&... inds)
    {
    return combiner(IndexSet(i1,std::forward<Inds>(inds)...));
    }

Index
combinedIndex(ITensor const& C);

//Construct diagonal ITensor with diagonal 
//elements set to 1.0
ITensor
delta(IndexSet const& is);

template<typename... Inds>
ITensor
delta(Index const& i1,
      Inds&&... inds)
    {
    return delta(IndexSet(i1,std::forward<Inds>(inds)...));
    }

//Construct diagonal ITensor,
//diagonal elements given by container C
//(Uses elements C.begin() up to C.end())
template<typename Container, 
         class = stdx::enable_if_t<stdx::containerOf<Real,Container>::value
                                || stdx::containerOf<Cplx,Container>::value> >
ITensor
diagITensor(Container const& C,
            IndexSet const& is);

template<typename Container, 
         typename... Inds,
         class = stdx::enable_if_t<stdx::containerOf<Real,Container>::value
                                || stdx::containerOf<Cplx,Container>::value> >
ITensor
diagITensor(Container const& C,
            Index const& i1,
            Inds&&... inds)
    {
    return diagITensor(C,IndexSet(i1,std::forward<Inds>(inds)...));
    }

//Construct ITensors with random elements
ITensor
randomITensor(IndexSet const& inds);
ITensor
randomITensorC(IndexSet const& inds);
ITensor
randomITensor(QN q, IndexSet const& inds);
ITensor
randomITensorC(QN q, IndexSet const& inds);

template <typename... Inds>
ITensor
randomITensor(Index const& i1,
              Inds&&... inds)
  {
  return randomITensor(IndexSet(i1,std::forward<Inds>(inds)...));
  }
template <typename... Inds>
ITensor
randomITensorC(Index const& i1,
               Inds&&... inds)
  {
  return randomITensorC(IndexSet(i1,std::forward<Inds>(inds)...));
  }
template <typename... Inds>
ITensor
randomITensor(QN q, Index const& i1, Inds&&... inds)
  {
  return randomITensor(q,IndexSet(i1,std::forward<Inds>(inds)...));
  }
template <typename... Inds>
ITensor
randomITensorC(QN q, Index const& i1, Inds&&... inds)
  {
  return randomITensorC(q,IndexSet(i1,std::forward<Inds>(inds)...));
  }

ITensor
matrixITensor(Matrix && M, IndexSet const& is);
ITensor
matrixITensor(Matrix const& M, IndexSet const& is);
ITensor
matrixITensor(CMatrix && M, IndexSet const& is);
ITensor
matrixITensor(CMatrix const& M, IndexSet const& is);

ITensor
matrixITensor(Matrix && M,
              Index const& i1, Index const& i2);
ITensor
matrixITensor(Matrix const& M,
              Index const& i1, Index const& i2);
ITensor
matrixITensor(CMatrix && M,
              Index const& i1, Index const& i2);
ITensor
matrixITensor(CMatrix const& M,
              Index const& i1, Index const& i2);

//
// QN ITensor related functions
//

QN
div(ITensor const& T);

//flux is an alias for div
QN
flux(ITensor const& T);

bool
hasQNs(ITensor const& T);

ITensor
toDense(ITensor T);

bool
isDense(ITensor const& T);

ITensor
removeQNs(ITensor T);

template<typename V>
TenRef<Range,V>
getBlock(ITensor & T, Block block_ind);

std::ostream& 
operator<<(std::ostream & s, ITensor const& T);

//
// Depecrated
//

void
randomize(ITensor & T, Args const& args = Args::global());

long
rank(ITensor const& T);

template<typename... Inds>
ITensor
reindex(ITensor const& cT, 
        Index o1, Index n1, 
        Inds... inds);

template<typename Container, 
         typename... Inds,
         class = stdx::enable_if_t<stdx::containerOf<Real,Container>::value
                                || stdx::containerOf<Cplx,Container>::value> >
ITensor
diagTensor(Container const& C,
           Index const& i1,
           Inds&&... inds);

template <typename... Inds>
ITensor
randomTensor(Index const& i1, Inds&&... inds);
template <typename... Inds>
ITensor
randomTensorC(Index const& i1, Inds&&... inds);
ITensor
randomTensor(IndexSet const& inds);

template <typename... Inds>
ITensor
randomTensor(QN q, Index const& i1, Inds&&... inds);
template <typename... Inds>
ITensor
randomTensorC(QN q, Index const& i1, Inds&&... inds);
ITensor
randomTensor(QN q, IndexSet const& inds, Args const& args = Args::global());

ITensor
matrixTensor(Matrix && M, Index const& i1, Index const& i2);
ITensor
matrixTensor(Matrix const& M, Index const& i1, Index const& i2);
ITensor
matrixTensor(CMatrix && M, Index const& i1, Index const& i2);
ITensor
matrixTensor(CMatrix const& M, Index const& i1, Index const& i2);

#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, ITensor const& I);
void
h5_read(h5::group parent, std::string const& name, ITensor & I);
#endif //ITENSOR_USE_HDF5

} //namespace itensor

#include "itensor_impl.h"


#endif
