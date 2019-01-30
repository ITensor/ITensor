//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_INTERFACE_IMPL_H_
#define __ITENSOR_ITENSOR_INTERFACE_IMPL_H_

#include "itensor/itdata/task_types.h"
#include "itensor/tensor/contract.h"
#include "itensor/iqindex.h"
//#include "itensor/util/print_macro.h"

//
// Template Method Implementations
//

namespace itensor {

QN
div(IQTensor const& T);

namespace detail {

void
allocReal(ITensor& T);

void
allocReal(IQTensor& T);

void
allocReal(ITensor& T, IntArray const& inds);

void
allocReal(IQTensor& T, IntArray const& inds);

void
allocCplx(ITensor& T);

void inline
allocCplx(IQTensor & T) { Error("allocCplx not defined for IQTensor"); }

} //namespace detail

template<typename IndexT>
template <typename... index_types>
ITensorT<IndexT>::
ITensorT(IndexT  const& i1,
         index_types const&... i2etc)
  : is_(i1,i2etc...)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

template<typename IndexT>
ITensorT<IndexT>::
ITensorT(std::vector<index_type> const& inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

template<typename IndexT>
template<size_t N> 
ITensorT<IndexT>::
ITensorT(std::array<index_type,N> const& inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

template<typename IndexT>
ITensorT<IndexT>::
ITensorT(std::initializer_list<index_type> inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

template<typename IndexT>
ITensorT<IndexT>::
ITensorT(indexset_type const& is)
  : is_(is)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

template<typename IndexT>
ITensorT<IndexT>::
ITensorT(indexset_type iset,
         storage_ptr&& pdat,
         LogNum const& scale)
    :
    is_(std::move(iset)),
    store_(std::move(pdat))
    { 
    IF_USESCALE(scale_ = scale;)
    }

template<typename IndexT>
template <class DataType>
ITensorT<IndexT>::
ITensorT(indexset_type iset,
         DataType&& dat,
         LogNum const& scale) :
    is_(std::move(iset)),
    store_(newITData<stdx::decay_t<DataType>>(std::move(dat)))
    {
    IF_USESCALE(scale_ = scale;)
    static_assert(std::is_rvalue_reference<decltype(std::forward<DataType>(dat))>::value,
                  "Error: cannot pass lvalues to ITensorT(...,DataType&& dat,...) constructor");
    }

template<typename IndexT>
ITensorT<IndexT>::
ITensorT(Cplx val) { Error("ITensorT(Cplx) not implemented"); }

template<typename IndexT>
ITensorT<IndexT>::
operator ITensor() const { Error("ITensorT->ITensor not implemented"); return *this; }

template<typename IndexT>
Cplx ITensorT<IndexT>::
cplx(std::vector<indexval_type> const& ivs) const
    {
    //using indexval_type = typename IndexT::indexval_type;

    if(!store()) Error("tensor storage unallocated");

    auto size = ivs.size();
    //constexpr size_t size = sizeof...(ivs)+1;
    //auto vals = std::array<indexval_type,size>{{static_cast<indexval_type>(iv1),static_cast<indexval_type>(ivs)...}};
    if(size != size_t(inds().r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ivs) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to real/cplx (expected %d, got %d)",inds().r(),size));
        }

    auto inds = IntArray(size);
    detail::permute_map(is_,ivs,inds,
                [](indexval_type const& iv) { return iv.val-1; });
    auto z = itensor::doTask(GetElt<IndexT>{is_,inds},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale().real0();
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in cplx(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in cplx(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

template<typename IndexT>
template<typename IV, typename... IVs>
auto ITensorT<IndexT>::
cplx(IV const& iv1, IVs&&... ivs) const
    -> stdx::if_compiles_return<Cplx,decltype(iv1.index),decltype(iv1.val)>
    {
    using indexval_type = typename IndexT::indexval_type;

    if(!store()) Error("tensor storage unallocated");

    constexpr size_t size = sizeof...(ivs)+1;
    auto vals = std::array<indexval_type,size>{{static_cast<indexval_type>(iv1),static_cast<indexval_type>(ivs)...}};
    if(size != size_t(inds().r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : vals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to real/cplx (expected %d, got %d)",inds().r(),size));
        }

    auto inds = IntArray(size);
    detail::permute_map(is_,vals,inds,
                [](indexval_type const& iv) { return iv.val-1; });
    auto z = itensor::doTask(GetElt<IndexT>{is_,inds},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale().real0(); 
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in cplx(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in cplx(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

template<typename IndexT>
template<typename Int>
auto ITensorT<IndexT>::
cplx(std::vector<Int> const& ints) const
    -> stdx::enable_if_t<std::is_integral<Int>::value,Cplx>
    {
    if(!store()) Error("tensor storage unallocated");

    auto size = ints.size();
    //constexpr size_t size = sizeof...(ivs)+1;
    //auto ints = std::array<Int,size>{{iv1,static_cast<int>(ivs)...}};
    if(size != size_t(inds().r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        print("Indices provided = ");
        for(auto i : ints) print(" ",i);
        println("\n---------------------------------------------");
        Error(format("Wrong number of ints passed to real/cplx (expected %d, got %d)",inds().r(),size));
        }

    auto inds = IntArray(size);
    for(auto i : range(size))
        inds[i] = ints[i]-1;
    auto z = itensor::doTask(GetElt<IndexT>{is_,inds},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale_.real0();
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in cplx(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in cplx(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

template<typename IndexT>
template<typename Int, typename... Ints>
auto ITensorT<IndexT>::
cplx(Int iv1, Ints... ivs) const
    -> stdx::enable_if_t<std::is_integral<Int>::value && stdx::and_<std::is_integral<Ints>...>::value,Cplx>
    {
    return this->cplx(std::vector<Int>{{iv1,static_cast<int>(ivs)...}});
    }

template<typename IndexT>
Cplx ITensorT<IndexT>::
cplx() const
    {
    if(inds().r() != 0)
        {
        Error(format("Wrong number of IndexVals passed to real/cplx (expected %d, got 0)",inds().r()));
        }
    constexpr size_t size = 0;
    auto inds = IntArray(size);
    auto z = itensor::doTask(GetElt<IndexT>{is_,inds},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale_.real0(); 
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in cplx(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in cplx(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

template<typename IndexT>
template <typename... IVals>
Real ITensorT<IndexT>::
real(IVals&&... ivs) const
    {
    auto z = cplx(std::forward<IVals>(ivs)...);
    if(fabs(z.imag()) > 1E-15 && fabs(z.imag()) > 1E-14*fabs(z.real()))
        {
        printfln("element = (%.5E,%.5E)",z.real(),z.imag());
        //Error("tensor is Complex valued, use .cplx(...) method");
        throw ITError("tensor is complex valued, use .cplx(...) method");
        }
    return z.real();
    }

namespace detail {

template<typename IndexValT, 
         typename Iter,
         typename IV>
auto
getVals(Iter it,
        Cplx & z,
        IV const& iv)
    -> stdx::if_compiles_return<void,decltype(iv.index),decltype(iv.val)>
    {
    static_assert(stdx::false_regardless_of<IndexValT>::value,
            "Last argument to .set method must be Real or Cplx scalar");
    }

template<typename IndexValT, typename Iter>
void
getVals(Iter it,
        Cplx & z,
        Cplx const& w)
    {
    z = w;
    }

template<typename IndexValT, 
         typename Iter, 
         typename IV, 
         typename... Rest>
auto
getVals(Iter it,
        Cplx & z,
        IV const& iv,
        Rest&&... rest)
    -> stdx::if_compiles_return<void,decltype(iv.index),decltype(iv.val)>
    {
    *it = static_cast<IndexValT>(iv);
    getVals<IndexValT>(++it,z,std::forward<Rest&&>(rest)...);
    }

template<typename IndexValT, typename Iter, typename... Rest>
bool
getVals(Iter it,
        Cplx & z,
        Cplx w,
        Rest&&... rest)
    {
    static_assert(stdx::false_regardless_of<Iter>::value,
            "New value passed to .set method must be last argument");
    return false;
    }

template<typename IntT,
         typename Iter,
         typename Arg>
auto
getInts(Iter it,
        Cplx & z,
        Arg const& iv)
    -> stdx::enable_if_t<not std::is_convertible<Arg,Cplx>::value,void>
    {
    static_assert(stdx::false_regardless_of<IntT>::value,
            "Last argument to .set method must be Real or Cplx scalar");
    }

template<typename IntT,
         typename Iter,
         typename Arg>
auto
getInts(Iter it,
        Cplx & z,
        Arg const& w)
    -> stdx::enable_if_t<std::is_convertible<Arg,Cplx>::value,void>
    {
    z = w;
    }

template<typename IntT,
         typename Iter,
         typename Int,
         typename... Rest>
auto
getInts(Iter it,
        Cplx & z,
        Int const& w,
        Rest&&... rest)
    -> stdx::enable_if_t<not std::is_integral<Int>::value,void>
    {
    static_assert(stdx::false_regardless_of<Iter>::value,
            "New value passed to .set method must be last argument");
    return false;
    }

template<typename IntT,
         typename Iter,
         typename Int,
         typename... Rest>
auto
getInts(Iter it,
        Cplx & z,
        Int w,
        Rest&&... rest)
    -> stdx::enable_if_t<std::is_integral<Int>::value,void>
    {
    *it = w-1;
    getInts<IntT>(++it,z,std::forward<Rest&&>(rest)...);
    }

} //namespace detail


template<typename IndexT>
template<typename IV, typename... VArgs>
auto ITensorT<IndexT>::
set(IV const& iv1, VArgs&&... vargs)
    -> stdx::if_compiles_return<void,decltype(iv1.index),decltype(iv1.val)>
    {
    static constexpr auto size = 1+(sizeof...(vargs)-1);
    std::array<indexval_type,size> vals;
    Cplx z;
    detail::getVals<indexval_type>(vals.begin(),z,iv1,std::forward<VArgs&&>(vargs)...);
    if(size != size_t(inds().r())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : vals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().r(),size));
        }
    auto inds = IntArray(is_.r(),0);
    detail::permute_map(is_,vals,inds,
                        [](indexval_type const& iv) { return iv.val-1; });
    //TODO: if !store_ and !is_real, call allocCplx instead
    //and move this line after check for is_real
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(z.imag()==0.0)
        {
        doTask(SetElt<Real,IndexT>{z.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx,IndexT>{z,is_,inds},store_);
        }
    }

template<typename IndexT>
template<typename Int, typename... VArgs>
auto ITensorT<IndexT>::
set(Int iv1, VArgs&&... vargs)
    -> stdx::enable_if_t<std::is_integral<Int>::value,void>
    {
    static constexpr auto size = 1+(sizeof...(vargs)-1);
    auto ints = IntArray(size,0);
    Cplx z;
    detail::getInts<Int>(ints.begin(),z,iv1,std::forward<VArgs&&>(vargs)...);
    if(size != size_t(inds().r())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        print("Indices provided =");
        for(auto& i : ints) print(" ",1+i);
        println();
        println("---------------------------------------------");
        Error(format("Wrong number of ints passed to set (expected %d, got %d)",
                     inds().r(),size));
        }
    //TODO: if !store_ and !is_real, call allocCplx instead
    //and move this line after check for is_real
    if(!store_) detail::allocReal(*this,ints);
    scaleTo(1.);
    if(z.imag()==0.0)
        {
        doTask(SetElt<Real,IndexT>{z.real(),is_,ints},store_);
        }
    else
        {
        doTask(SetElt<Cplx,IndexT>{z,is_,ints},store_);
        }
    }

template<typename IndexT>
void ITensorT<IndexT>::
set(Cplx val)
    {
    if(0 != size_t(inds().r())) 
        {
        Error(format("Wrong number of IndexVals passed to set (expected %d, got 0)",
                     inds().r()));
        }
    auto inds = IntArray(0,1);
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.)
        {
        doTask(SetElt<Real,IndexT>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx,IndexT>{val,is_,inds},store_);
        }
    }

template<typename IndexT>
void ITensorT<IndexT>::
set(std::vector<typename IndexT::indexval_type> const& ivals,
    Cplx val)
    {
    auto size = ivals.size();
    if(size != size_t(inds().r())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ivals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().r(),size));
        }
    auto inds = IntArray(is_.r(),0);
    detail::permute_map(is_,ivals,inds,
                        [](indexval_type const& iv) { return iv.val-1; });
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.0)
        {
        doTask(SetElt<Real,IndexT>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx,IndexT>{val,is_,inds},store_);
        }
    }

template<typename IndexT>
void ITensorT<IndexT>::
set(std::vector<int> const& ints,
    Cplx val)
    {
    auto size = ints.size();
    if(size != size_t(inds().r())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ints) println(iv);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().r(),size));
        }
    auto inds = IntArray(is_.r(),0);
    for(auto i : range(size))
        inds[i] = ints[i]-1;
    //TODO: if !store_ and !is_real, call allocCplx instead
    //detail::permute_map(is_,ivals,inds,
    //                    [](indexval_type const& iv) { return iv.val-1; });
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.0)
        {
        doTask(SetElt<Real,IndexT>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx,IndexT>{val,is_,inds},store_);
        }
    }

template<typename IndexT>
template <typename Func>
ITensorT<IndexT>& ITensorT<IndexT>::
generate(Func&& f)
    {
    if(not this->store())
        {
        using RetType = decltype(f());
        if(std::is_same<RetType,Real>::value)
            {
            detail::allocReal(*this); 
            }
        else if(std::is_same<RetType,Cplx>::value)
            {
            detail::allocCplx(*this); 
            }
        else
            {
            Error("generate: generator function must return Real or Cplx scalar value");
            }
        }
    scaleTo(1);
    doTask(GenerateIT<decltype(f)>{std::forward<Func>(f)},store_);
    return *this;
    }

template<typename IndexT>
template <typename Func>
ITensorT<IndexT>& ITensorT<IndexT>::
apply(Func&& f)
    {
    scaleTo(1);
    doTask(ApplyIT<decltype(f)>{std::forward<Func>(f)},store_);
    return *this;
    }

template<typename IndexT>
template <typename Func>
const ITensorT<IndexT>& ITensorT<IndexT>::
visit(Func&& f) const
    {
    doTask(VisitIT<decltype(f)>{std::forward<Func>(f),LogNum{scale().real0()}},store_);
    return *this;
    }

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
conj()
    {
    doTask(Conj{},store_);
    return *this;
    }

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
dag()
    {
    Error("dag not implemented");
    return *this;
    }

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
takeReal()
    {
    doTask(TakeReal{},store_);
    return *this;
    }

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
takeImag()
    {
    doTask(TakeImag{},store_);
    return *this;
    }

#ifdef USESCALE
template<typename IndexT> 
void ITensorT<IndexT>::
scaleTo(scale_type const& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    doTask(Mult<Real>{scale_.real0()},store_);
    scale_ = newscale;
    }

template<typename IndexT> 
void ITensorT<IndexT>::
scaleTo(Real newscale) { scaleTo(LogNum{newscale}); }
#endif

template<typename IndexT>
void ITensorT<IndexT>::
swap(ITensorT & other)
    {
    is_.swap(other.is_);
    store_.swap(other.store_);
    IF_USESCALE(scale_.swap(other.scale_);)
    }

template <typename IVal, typename... IVals>
ITensorT<typename std::common_type<IVal,IVals...>::type::index_type>
setElt(IVal const& iv1, 
       IVals const&... rest)
    {
    //using index_type = typename IVal::index_type;
    using index_type = typename std::common_type<IVal,IVals...>::type::index_type;
    const constexpr auto size = 1+sizeof...(rest);
    auto ivs = stdx::make_array(iv1,rest...);
    //TODO: try directly making inds as iv1.index,(rest.index)...
    auto inds = std::array<index_type,size>{};
    for(size_t j = 0; j < size; ++j) inds[j] = ivs[j].index;
    auto D = ITensorT<index_type>{IndexSetT<index_type>(inds)};
    D.set(iv1,rest...,1.);
    return D;
    }

#ifndef USESCALE

template <typename IndexT>
ITensorT<IndexT> ITensorT<IndexT>::
operator-() const
    { 
    auto res = *this;
    doTask(Mult<Real>(-1.),res.store());
    return res;
    }

#else

template <typename IndexT>
ITensorT<IndexT> ITensorT<IndexT>::
operator-() const
    { 
    auto res = *this;
    res.scale_.negate();
    return res;
    }

#endif

template<typename I>
ITensorT<I> 
operator*(ITensorT<I> A, ITensorT<I> const& B) { A *= B; return A; }
template<typename I>
ITensorT<I>
operator*(ITensorT<I> const& A, ITensorT<I>&& B) { B *= A; return B; }
template<typename I>
ITensorT<I> 
operator*(ITensorT<I> T, Real fac) { T *= fac; return T; }
template<typename I>
ITensorT<I> 
operator*(Real fac, ITensorT<I> T) { T *= fac; return T; }
template<typename I>
ITensorT<I> 
operator*(ITensorT<I> T, Complex fac) { T *= fac; return T; }
template<typename I>
ITensorT<I> 
operator*(Complex fac, ITensorT<I> T) { T *= fac; return T; }
template<typename I>
ITensorT<I> 
operator/(ITensorT<I> T, Real fac) { T /= fac; return T; }
template<typename I>
ITensorT<I> 
operator/(ITensorT<I> T, Complex fac) { T /= fac; return T; }

template<typename I>
ITensorT<I> 
operator+(ITensorT<I> A, const ITensorT<I>& B) { A += B; return A; }
template<typename I>
ITensorT<I> 
operator+(ITensorT<I> const& A, ITensorT<I>&& B) { B += A; return B; }

template<typename I>
ITensorT<I> 
operator-(ITensorT<I> A, const ITensorT<I>& B) { A -= B; return A; }
template<typename I>
ITensorT<I> 
operator-(ITensorT<I> const& A, ITensorT<I>&& B) { B -= A; B *= -1; return B; }

template<typename I>
ITensorT<I> 
operator/(ITensorT<I> A, ITensorT<I> const& B) { A /= B; return A; }
template<typename I>
ITensorT<I>
operator/(ITensorT<I> const& A, ITensorT<I> && B) { B /= A; return B; }


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
primeLevel(ITensorT<IndexT> A, 
           VarArgs&&... vargs)
    {
    A.primeLevel(std::forward<VarArgs>(vargs)...);
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

template<typename IndexT, typename... VarArgs>
ITensorT<IndexT>
sim(ITensorT<IndexT> A, 
    VarArgs&&... vargs)
    {
    A.sim(std::forward<VarArgs>(vargs)...);
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

template<typename IndexT,
         typename Cond>
IndexT
findindex(ITensorT<IndexT> const& T, 
          Cond && cond)
    {
    for(auto& i : T.inds())
        if(cond(i)) return i;
    return IndexT{};
    }

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
template<typename IndexT> 
IndexT
commonIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t);


//Find index of tensor A (of optional type t) 
//which is NOT shared by tensor B
template<typename IndexT> 
IndexT
uniqueIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t);

template<typename IndexT, typename... Tensors> 
IndexT
uniqueIndex(ITensorT<IndexT> const& A, 
            ITensorT<IndexT> const& T1,
            ITensorT<IndexT> const& T2,
            Tensors const&... Tens)
    {
    auto Ts = stdx::make_array(T1,T2,Tens...);
    for(auto& I : A.inds())
        {
        bool found = false;
        for(auto& T : Ts) if(hasindex(T,I))
            {
            found = true;
            break;
            }
        if(!found) return I;
        }
    return IndexT();
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
          IndexType type);

//Apply x = f(x) for each element x of T
//and return the resulting tensor
template<typename I, typename F>
ITensorT<I>
apply(ITensorT<I> T, F&& f)
    {
    T.apply(std::forward<F>(f));
    return T;
    }

template<typename I>
bool
isComplex(ITensorT<I> const& T)
    {
    return doTask(CheckComplex{},T.store());
    }

template<typename I>
bool
isReal(ITensorT<I> const& T)
    {
    return not isComplex(T);
    }

template<typename I>
long
rank(ITensorT<I> const& T) { return rank(T.inds()); }

//return number of indices of T
//(same as rank)
template<typename I>
long
ord(ITensorT<I> const& T) { return rank(T.inds()); }

template<typename I>
Real
norm(ITensorT<I> const& T);

template<typename I>
void
randomize(ITensorT<I> & T, Args const& args);

template<typename I>
ITensorT<I>
random(ITensorT<I> T, const Args& args)
    {
    randomize(T,args);
    return T;
    }

template<typename I>
ITensorT<I>
conj(ITensorT<I> T)
    {
    T.conj();
    return T;
    }

template<typename I>
ITensorT<I>
dag(ITensorT<I> T)
    {
    T.dag();
    return T;
    }

template<typename I>
Real
sumels(const ITensorT<I>& t)
    {
    auto z = sumelsC(t);
    if(z.imag() != 0) Error("ITensor has non-zero imaginary part, use sumelsC");
    return z.real();
    }

template<typename I>
Cplx
sumelsC(ITensorT<I> const& t)
    {
    auto z = doTask(SumEls<I>{t.inds()},t.store());
#ifndef USESCALE
    return z;
#else
    return t.scale().real0()*z;
#endif
    }

template<typename IndexT, typename... Inds>
ITensorT<IndexT>
reindex(ITensorT<IndexT> const& cT, 
        IndexT o1, IndexT n1, 
        Inds... inds) 
    {
    constexpr size_t size = 2+sizeof...(inds);
    auto ipairs = std::array<IndexT,size>{{o1,n1,static_cast<IndexT>(inds)...}};

    auto T = cT;
    auto is = T.inds();

    for(auto j : range(is))
        {
        for(size_t oi = 0, ni = 1; ni <= size; oi += 2, ni += 2)
            {
            if(is[j].noprimeEquals(ipairs[oi]))
                {
                if(is[j].m() != ipairs[ni].m())
                    {
                    printfln("Old m = %d",is[j].m());
                    printfln("New m would be = %d",ipairs[ni].m());
                    throw ITError("Mismatch of index dimension in reindex");
                    }
                auto plev = is[j].primeLevel();
                auto arrow_dir = is[j].dir();
                is[j] = noprime(ipairs[ni]);
                is[j].primeLevel(plev);
                is[j].dir(arrow_dir);
                break;
                }
            }
        }
    auto nT = ITensorT<IndexT>(is,std::move(T.store()),T.scale());
    return nT;
    }


namespace detail {


// TODO: check for incorrect inputs

template<typename IndexT,
         typename Iter>
void
getDotInds(Iter it,
           std::string const& dots)
    {
    if(dots != "...")
        Error(format("Wrong string passed to order (expected '...', got '%s')",dots));
    }

template<typename IndexT,
         typename Iter,
         typename... Rest>
void
getDotInds(Iter it,
           IndexT const& ind,
           Rest const&... rest)
    {
    *it = ind;
    getDotInds<IndexT>(++it,std::forward<Rest const&>(rest)...);
    }

template <typename IndexT>
IndexSetT<IndexT>
moveToFront(IndexSetT<IndexT> const& isf, IndexSetT<IndexT> const& is)
    {
    auto rf = isf.r();
    auto r = is.r();

    if(rf >= r)
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n",isf," '...'\n");
        println("---------------------------------------------");
        Error(format("Wrong number of indices passed to order (expected < %d, got %d)",r,rf));
        }

    auto iso = IndexSetT<IndexT>(r);

    auto i = 0;
    for(auto& I : isf) 
        {
        if(!hasindex(is,I))
            {
            println("---------------------------------------------");
            println("Tensor indices = \n",is,"\n");
            println("---------------------------------------------");
            println("Indices provided = \n",isf," '...'\n");
            println("---------------------------------------------");
            Error(format("Bad index passed to order"));
            }
        iso[i] = I;
        i++;
        }

    auto j = rf;
    for(auto& J : is)
        {
        if(!hasindex(isf,J))
            {
            iso[j] = J;
            j++;
            }
        }

    return iso;
    }

template <typename IndexT>
IndexSetT<IndexT>
moveToBack(IndexSetT<IndexT> const& isb, IndexSetT<IndexT> const& is)
    {
    auto rb = isb.r();
    auto r = is.r();

    if(rb >= r)
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n'...' ",isb,"\n");
        println("---------------------------------------------");
        Error(format("Wrong number of indices passed to order (expected < %d, got %d)",r,rb));
        }

    auto iso = IndexSetT<IndexT>(r);

    auto i = r-rb;
    for(auto& I : isb) 
        {
        if(!hasindex(is,I))
            {
            println("---------------------------------------------");
            println("Tensor indices = \n",is,"\n");
            println("---------------------------------------------");
            println("Indices provided = \n'...' ",isb,"\n");
            println("---------------------------------------------");
            Error(format("Bad index passed to order"));
            }
        iso[i] = I;
        i++;
        }

    auto j = 0;
    for(auto& J : is)
        {
        if(!hasindex(isb,J))
            {
            iso[j] = J;
            j++;
            }
        }

    return iso;
    }

} //namespace detail

//Version of order accepting syntax: T.order(i,j,"...")
template <typename IndexT>
template <typename... Indxs>
auto ITensorT<IndexT>::
order(IndexT const& ind1, Indxs const&... inds)
    -> stdx::enable_if_t<not stdx::and_<std::is_same<IndexT, Indxs>...>::value,ITensorT<IndexT>&>
    {
    static constexpr auto size = 1+(sizeof...(inds)-1);
    auto isf = IndexSetT<IndexT>(size);
    detail::getDotInds<IndexT>(isf.begin(),ind1,std::forward<Indxs const&>(inds)...);
    order(detail::moveToFront(isf,this->inds()));
    return *this;
    }

//Version of order accepting syntax: T.order(i,j,k)
template <typename IndexT>
template <typename... Indxs>
auto ITensorT<IndexT>::
order(IndexT const& ind1, Indxs const&... inds)
    -> stdx::enable_if_t<stdx::and_<std::is_same<IndexT, Indxs>...>::value,ITensorT<IndexT>&>
    {
    order(IndexSetT<IndexT>(ind1, inds...));
    return *this;
    }

//Version of order accepting syntax: T.order("...",j,k)
template <typename IndexT>
template <typename... Indxs>
ITensorT<IndexT>& ITensorT<IndexT>::
order(std::string const& dots, Indxs const&... inds)
    {
    if(dots != "...")
        Error(format("Wrong string passed to order (expected '...', got '%s')",dots));
    order(detail::moveToBack(IndexSetT<IndexT>(inds...),this->inds()));
    return *this;
    }

//order function which returns a new ITensor 
template <typename IndexT, typename... Indxs>
ITensorT<IndexT>
order(ITensorT<IndexT> A, Indxs const&... inds)
    {
    A.order(std::forward<Indxs const&>(inds)...);
    return A;
    }

template<typename T, typename... CtrArgs>
ITensor::storage_ptr
readType(std::istream& s, CtrArgs&&... args)
    {
    T t(std::forward<CtrArgs>(args)...);
    read(s,t);
    return newITData<T>(std::move(t));
    }

template<typename I>
void ITensorT<I>::
read(std::istream& s)
    {
    itensor::read(s,is_);
    LogNum scale;
    itensor::read(s,scale);
    IF_USESCALE(scale_ = scale;)
    auto type = StorageType::Null;
    itensor::read(s,type);
    if(type==StorageType::Null) { /*intentionally left blank*/  }
    else if(type==StorageType::DenseReal) { store_ = readType<DenseReal>(s); }
    else if(type==StorageType::DenseCplx) { store_ = readType<DenseCplx>(s); }
    else if(type==StorageType::Combiner) { store_ = readType<Combiner>(s); }
    else if(type==StorageType::DiagReal) { store_ = readType<Diag<Real>>(s); }
    else if(type==StorageType::DiagCplx) { store_ = readType<Diag<Cplx>>(s); }
    else if(type==StorageType::QDenseReal) { store_ = readType<QDense<Real>>(s); }
    else if(type==StorageType::QDenseCplx) { store_ = readType<QDense<Cplx>>(s); }
    else if(type==StorageType::QDiagReal) { store_ = readType<QDiag<Real>>(s); }
    else if(type==StorageType::QDiagCplx) { store_ = readType<QDiag<Cplx>>(s); }
    else if(type==StorageType::QCombiner) { store_ = readType<QCombiner>(s); }
    else if(type==StorageType::ScalarReal) { store_ = readType<ScalarReal>(s); }
    else if(type==StorageType::ScalarCplx) { store_ = readType<ScalarCplx>(s); }
    else
        {
        Error("Unrecognized type when reading tensor from istream");
        }
    }

struct Write
    {
    std::ostream& s;

    Write(std::ostream& s_) : s(s_) { }
    };

inline const char*
typeNameOf(Write const&) { return "Write"; }

template<typename D>
auto
doTask(Write & W, D const& d)
    -> stdx::if_compiles_return<void,decltype(itensor::write(W.s,d))>
    {
    write(W.s,d);
    }

//template<typename I>
//void
//write(std::ostream& s, ITensorT<I> const& T);


} // namespace itensor


#endif
