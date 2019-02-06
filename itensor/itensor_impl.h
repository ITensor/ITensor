//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_INTERFACE_IMPL_H_
#define __ITENSOR_ITENSOR_INTERFACE_IMPL_H_

#include "itensor/itdata/task_types.h"
#include "itensor/tensor/contract.h"
//#include "itensor/util/print_macro.h"

//
// Template Method Implementations
//

namespace itensor {


namespace detail {

void
allocReal(ITensor& T);

void
allocReal(ITensor& T, IntArray const& inds);

void
allocCplx(ITensor& T);

} //namespace detail

template <typename... index_types>
ITensor::
ITensor(Index  const& i1,
        index_types const&... i2etc)
  : is_(i1,i2etc...)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }


template<size_t N> 
ITensor::
ITensor(std::array<Index,N> const& inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }




template <class DataType>
ITensor::
ITensor(IndexSet iset,
        DataType&& dat,
        LogNum const& scale) :
    is_(std::move(iset)),
    store_(newITData<stdx::decay_t<DataType>>(std::move(dat)))
    {
    IF_USESCALE(scale_ = scale;)
    static_assert(std::is_rvalue_reference<decltype(std::forward<DataType>(dat))>::value,
                  "Error: cannot pass lvalues to ITensor(...,DataType&& dat,...) constructor");
    }

template<typename IV, typename... IVs>
auto ITensor::
cplx(IV const& iv1, IVs&&... ivs) const
    -> stdx::if_compiles_return<Cplx,decltype(iv1.index),decltype(iv1.val)>
    {
    constexpr size_t size = sizeof...(ivs)+1;
    auto vals = std::array<IndexVal,size>{{static_cast<IndexVal>(iv1),
                                           static_cast<IndexVal>(ivs)...}};
    if(size != size_t(inds().r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : vals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to real/cplx (expected %d, got %d)",
                     inds().r(),size));
        }

    auto inds = IntArray(size);
    detail::permute_map(is_,vals,inds,
                [](IndexVal const& iv) { return iv.val-1; });
    auto z = itensor::doTask(GetElt{is_,inds},store_);
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


template<typename Int>
auto ITensor::
cplx(std::vector<Int> const& ints) const
    -> stdx::enable_if_t<std::is_integral<Int>::value,Cplx>
    {
    if(!store()) Error("tensor storage unallocated");

    auto size = ints.size();
    if(size != size_t(inds().r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        print("Indices provided = ");
        for(auto i : ints) print(" ",i);
        println("\n---------------------------------------------");
        Error(format("Wrong number of ints passed to real/cplx (expected %d, got %d)",
                     inds().r(),size));
        }

    auto inds = IntArray(size);
    for(auto i : range(size))
        inds[i] = ints[i]-1;
    auto z = itensor::doTask(GetElt{is_,inds},store_);
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

template<typename Int, typename... Ints>
auto ITensor::
cplx(Int iv1, Ints... ivs) const
    -> stdx::enable_if_t<std::is_integral<Int>::value 
                     && stdx::and_<std::is_integral<Ints>...>::value,Cplx>
    {
    return this->cplx(std::vector<Int>{{iv1,static_cast<int>(ivs)...}});
    }

template <typename... IVals>
Real ITensor::
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


template<typename IV, typename... VArgs>
auto ITensor::
set(IV const& iv1, VArgs&&... vargs)
    -> stdx::if_compiles_return<void,decltype(iv1.index),decltype(iv1.val)>
    {
    static constexpr auto size = 1+(sizeof...(vargs)-1);
    std::array<IndexVal,size> vals;
    Cplx z;
    detail::getVals<IndexVal>(vals.begin(),z,iv1,std::forward<VArgs&&>(vargs)...);
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
                        [](IndexVal const& iv) { return iv.val-1; });
    //TODO: if !store_ and !is_real, call allocCplx instead
    //and move this line after check for is_real
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(z.imag()==0.0)
        {
        doTask(SetElt<Real>{z.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{z,is_,inds},store_);
        }
    }

template<typename Int, typename... VArgs>
auto ITensor::
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
        doTask(SetElt<Real>{z.real(),is_,ints},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{z,is_,ints},store_);
        }
    }



template <typename Func>
ITensor& ITensor::
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

template <typename Func>
ITensor& ITensor::
apply(Func&& f)
    {
    scaleTo(1);
    doTask(ApplyIT<decltype(f)>{std::forward<Func>(f)},store_);
    return *this;
    }

template <typename Func>
const ITensor& ITensor::
visit(Func&& f) const
    {
    doTask(VisitIT<decltype(f)>{std::forward<Func>(f),LogNum{scale().real0()}},store_);
    return *this;
    }


#ifdef USESCALE
void inline ITensor::
scaleTo(scale_type const& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    doTask(Mult<Real>{scale_.real0()},store_);
    scale_ = newscale;
    }

void inline ITensor::
scaleTo(Real newscale) { scaleTo(LogNum{newscale}); }
#endif

template <typename IVal, typename... IVals>
ITensor
setElt(IVal const& iv1, 
       IVals const&... rest)
    {
    using index_type = typename std::common_type<IVal,IVals...>::type::index_type;
    const constexpr auto size = 1+sizeof...(rest);
    auto ivs = stdx::make_array(iv1,rest...);
    //TODO: try directly making inds as iv1.index,(rest.index)...
    auto inds = std::array<index_type,size>{};
    for(size_t j = 0; j < size; ++j) inds[j] = ivs[j].index;
    auto D = ITensor{IndexSet(inds)};
    D.set(iv1,rest...,1.);
    return D;
    }

#ifndef USESCALE

ITensor inline ITensor::
operator-() const
    { 
    auto res = *this;
    doTask(Mult<Real>(-1.),res.store());
    return res;
    }

#else

ITensor inline ITensor::
operator-() const
    { 
    auto res = *this;
    res.scale_.negate();
    return res;
    }

#endif

ITensor inline
operator*(ITensor A, ITensor const& B) { A *= B; return A; }

ITensor inline
operator*(ITensor const& A, ITensor&& B) { B *= A; return B; }
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
operator+(ITensor A, ITensor const& B) { A += B; return A; }
ITensor inline
operator+(ITensor const& A, ITensor&& B) { B += A; return B; }
ITensor inline
operator-(ITensor A, ITensor const& B) { A -= B; return A; }
ITensor inline
operator-(ITensor const& A, ITensor&& B) { B -= A; B *= -1; return B; }
ITensor inline
operator/(ITensor A, ITensor const& B) { A /= B; return A; }
ITensor inline
operator/(ITensor const& A, ITensor && B) { B /= A; return B; }


template<typename... VarArgs>
ITensor
prime(ITensor A, 
      VarArgs&&... vargs)
    {
    A.prime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
setPrime(ITensor A,
         VarArgs&&... vargs)
    {
    A.setPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
noPrime(ITensor A, 
        VarArgs&&... vargs)
    {
    A.noPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
mapPrime(ITensor A, 
         VarArgs&&... vargs)
    {
    A.mapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
swapPrime(ITensor A, 
          VarArgs&&... vargs)
    {
    A.swapPrime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
replaceTags(ITensor A,
            VarArgs&&... vargs)
    {
    A.replaceTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
setTags(ITensor A,
        VarArgs&&... vargs)
    {
    A.setTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
addTags(ITensor A,
        VarArgs&&... vargs)
    {
    A.addTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
removeTags(ITensor A,
           VarArgs&&... vargs)
    {
    A.removeTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
swapTags(ITensor A,
         VarArgs&&... vargs)
    {
    A.swapTags(std::forward<VarArgs>(vargs)...);
    return A;
    }

//TODO: bring this back?
//template<typename Cond>
//Index
//findIndex(ITensor const& T, 
//          Cond && cond)
//    {
//    for(auto& i : T.inds()) if(cond(i)) return i;
//    return Index{};
//    }

Index inline
findIndex(ITensor const& T,
          TagSet const& tsmatch, 
          int plmatch)
    {
    return findIndex(T.inds(),tsmatch,plmatch);
    }

Index inline
findIndexExact(ITensor const& T,
               TagSet const& tsmatch, 
               int plmatch)
    {
    return findIndex(T.inds(),tsmatch,plmatch);
    }

//Apply x = f(x) for each element x of T
//and return the resulting tensor
template<typename F>
ITensor
apply(ITensor T, F&& f)
    {
    T.apply(std::forward<F>(f));
    return T;
    }


long inline
rank(ITensor const& T) { return rank(T.inds()); }

//return number of indices of T
//(same as rank)
long inline
ord(ITensor const& T) { return rank(T.inds()); }

Real
norm(ITensor const& T);

void
randomize(ITensor & T, Args const& args);

template <typename... Inds>
ITensor
randomITensor(Index const& i1, Inds&&... inds)
    {
    return random(ITensor(i1,std::forward<Inds>(inds)...));
    }
template <typename... Inds>
ITensor
randomITensorC(Index const& i1, Inds&&... inds)
    {
    return random(ITensor(i1,std::forward<Inds>(inds)...),{"Complex",true});
    }

template <typename... Inds>
ITensor
randomITensor(QN q, Index const& i1, Inds&&... inds)
    {
    auto is = IndexSet(i1,std::forward<Inds>(inds)...);
    return randomITensor(q,std::move(is));
    }
template <typename... Inds>
ITensor
randomITensorC(QN q, Index const& i1, Inds&&... inds)
    {
    auto is = IndexSet(i1,std::forward<Inds>(inds)...);
    return randomITensor(q,std::move(is),{"Complex=",true});
    }

//Deprecated
template <typename... Inds>
ITensor
randomTensor(Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensor(Index,...) is deprecated in favor of randomITensor(Index,...)");
    return randomITensor(i1,inds...);
    }
template <typename... Inds>
ITensor
randomTensorC(Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensorC(Index,...) is deprecated in favor of randomITensorC(Index,...)");
    return randomITensorC(i1,inds...);
    }

//Deprecated
template <typename... Inds>
ITensor
randomTensor(QN q, Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensor(QN,Index,...) is deprecated in favor of randomITensor(QN,Index,...)");
    return randomITensor(q,i1,inds...);
    }
template <typename... Inds>
ITensor
randomTensorC(QN q, Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensorC(QN,Index,...) is deprecated in favor of randomITensorC(QN,Index,...)");
    return randomITensorC(q,i1,inds...);
    }


ITensor inline
conj(ITensor T)
    {
    T.conj();
    return T;
    }

ITensor inline
dag(ITensor T)
    {
    T.dag();
    return T;
    }


template<typename... Inds>
ITensor
replaceInds(ITensor const& cT, 
            Index o1, Index n1, 
            Inds... inds) 
    {
    constexpr size_t size = 2+sizeof...(inds);
    auto ipairs = std::array<Index,size>{{o1,n1,static_cast<Index>(inds)...}};

    auto T = cT;
    auto is = T.inds();
    
    //This is a random prime level increase to 
    //prevent clashing indices if there are prime level
    //swaps
    auto tempLevel = 43218;
    for(auto J : is)
        {
        for(size_t oi = 0, ni = 1; ni <= size; oi += 2, ni += 2)
            {
            if(J==ipairs[oi])
                {
                if(J.m() != ipairs[ni].m())
                    {
                    printfln("Old m = %d",J.m());
                    printfln("New m would be = %d",ipairs[ni].m());
                    throw ITError("Mismatch of index dimension in reindex");
                    }
                T *= delta(J,prime(ipairs[ni],tempLevel));
                break;
                }
            }
        }

    //Bring the prime levels back down to the original
    //desired ones
    for(size_t ni = 1; ni <= size; ni += 2) T.prime(-tempLevel,prime(ipairs[ni],tempLevel));

    return T;
    }

template<typename... Inds>
ITensor
reindex(ITensor const& cT, 
        Index o1, Index n1, 
        Inds... inds) 
    {
    Global::warnDeprecated("reindex(ITensor,Index,Index,...) is deprecated in favor of replaceInds(ITensor,Index,Index,...)");
    constexpr size_t size = 2+sizeof...(inds);
    auto ipairs = std::array<Index,size>{{o1,n1,static_cast<Index>(inds)...}};

    auto T = cT;
    auto is = T.inds();

    for(auto j : range(is))
        {
        for(size_t oi = 0, ni = 1; ni <= size; oi += 2, ni += 2)
            {
            if(equalsIgnorePrime(is[j],ipairs[oi]))
                {
                if(is[j].m() != ipairs[ni].m())
                    {
                    printfln("Old m = %d",is[j].m());
                    printfln("New m would be = %d",ipairs[ni].m());
                    throw ITError("Mismatch of index dimension in reindex");
                    }
                auto plev = primeLevel(is[j]);
                auto arrow_dir = is[j].dir();
                is[j] = noPrime(ipairs[ni]);
                is[j].setPrime(plev);
                is[j].dir(arrow_dir);
                break;
                }
            }
        }
    auto nT = ITensor(is,std::move(T.store()),T.scale());
    return nT;
    }


namespace detail {


// TODO: check for incorrect inputs

template<typename Iter>
void
getDotInds(Iter it,
           std::string const& dots)
    {
    if(dots != "...")
        Error(format("Wrong string passed to order (expected '...', got '%s')",dots));
    }

template<typename Iter,
         typename... Rest>
void
getDotInds(Iter it,
           Index const& ind,
           Rest const&... rest)
    {
    *it = ind;
    getDotInds(++it,std::forward<Rest const&>(rest)...);
    }

IndexSet
moveToFront(IndexSet const& isf, IndexSet const& is);

IndexSet 
moveToBack(IndexSet const& isb, IndexSet const& is);


} //namespace detail

//Version of permute accepting syntax: T.permute(i,j,"...")
template <typename... Indxs>
auto ITensor::
permute(Index const& ind1, Indxs const&... inds)
        -> stdx::enable_if_t<not stdx::and_<std::is_same<Index, Indxs>...>::value,ITensor&>
    {
    static constexpr auto size = 1+(sizeof...(inds)-1);
    auto isf = IndexSet(size);
    detail::getDotInds(isf.begin(),ind1,std::forward<Indxs const&>(inds)...);
    permute(detail::moveToFront(isf,this->inds()));
    return *this;
    }

//Version of order accepting syntax: T.permute(i,j,k)
template <typename... Indxs>
auto ITensor::
permute(Index const& ind1, Indxs const&... inds)
        -> stdx::enable_if_t<stdx::and_<std::is_same<Index, Indxs>...>::value,ITensor&>
    {
    permute(IndexSet(ind1, inds...));
    return *this;
    }

//Version of order accepting syntax: T.permute("...",j,k)
template <typename... Indxs>
ITensor& ITensor::
permute(std::string const& dots, Indxs const&... inds)
    {
    if(dots != "...")
        Error(format("Wrong string passed to order (expected '...', got '%s')",dots));
    permute(detail::moveToBack(IndexSet(inds...),this->inds()));
    return *this;
    }

//order function which returns a new ITensor 
template<typename... Indxs>
ITensor
permute(ITensor A, Indxs const&... inds)
    {
    A.permute(std::forward<Indxs const&>(inds)...);
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
//write(std::ostream& s, ITensor<I> const& T);


template<typename... Inds>
ITensor
delta(Index const& i1,
      Inds const&... inds)
    { 
    auto is = IndexSet(i1,inds...);
    if(hasQNs(is))
        {
        return ITensor(std::move(is),QDiagReal(is,1.));
        }
    auto len = minM(is);
    return ITensor(std::move(is),DiagReal(len,1.));
    }

template<typename Container, typename... Inds, class>
ITensor
diagITensor(Container const& C, 
            Index const& i1,
            Inds &&... inds)
    { 
    auto is = IndexSet(i1,std::forward<Inds>(inds)...);
#ifdef DEBUG
    using size_type = decltype(C.size());
    //Compute min of all index dimensions
    auto minm = i1.m();
    for(const auto& ind : is)
        if(ind.m() < minm) minm = ind.m();
    if(C.size() != size_type(minm))
        {
        println("minm = ",minm);
        println("C.size() = ",C.size());
        Error("Wrong size of data in diagonal ITensor constructor");
        }
#endif
    using value_type = typename Container::value_type;
    return ITensor(std::move(is),Diag<value_type>(C.begin(),C.end()));
    }

//Deprecated
template<typename Container, typename... Inds, class>
ITensor
diagTensor(Container const& C, 
           Index const& i1,
           Inds &&... inds)
    { 
    Global::warnDeprecated("diagTensor(Container,Index,...) is deprecated in favor of diagITensor(Container,Index,...)");
    return diagITensor(C,i1,inds...);
    }

bool inline
hasQNs(ITensor const& T) { return hasQNs(T.inds()); }

template<typename V>
TenRef<Range,V>
getBlock(ITensor & T,
         IntArray block_ind)
    {
    if(block_ind.size() != size_t(T.r())) Error("Mismatched number of indices and ITensor rank");
    if(not T.store())
        {
        QN q;
        for(auto n : range(block_ind))
            {
            auto& I = T.inds()[n];
            q += I.qn(block_ind[n])*I.dir();
            }
        T = ITensor(T.inds(),QDense<V>(T.inds(),q));
        }
    //Interface is 1-indexed; switch to 0-indexed
    for(auto& i : block_ind) { i -= 1; }
    auto G = GetBlock<V>(T.inds(),block_ind);
    return doTask(G,T.store());
    }

template<typename... Tensors> 
Index
uniqueIndex(ITensor const& A, 
            ITensor const& T1,
            ITensor const& T2,
            Tensors const&... Tens)
    {
    auto Ts = stdx::make_array(T1,T2,Tens...);
    for(auto& I : A.inds())
        {
        bool found = false;
        for(auto& T : Ts) if(hasIndex(T,I))
            {
            found = true;
            break;
            }
        if(!found) return I;
        }
    return Index();
    }

} // namespace itensor


#endif
