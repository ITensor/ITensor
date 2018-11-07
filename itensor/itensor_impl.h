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
ITensor(std::array<index_type,N> const& inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }




template <class DataType>
ITensor::
ITensor(indexset_type iset,
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

    if(!store()) Error("tensor storage unallocated");

    constexpr size_t size = sizeof...(ivs)+1;
    auto vals = std::array<IndexVal,size>{{static_cast<IndexVal>(iv1),static_cast<IndexVal>(ivs)...}};
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

template<typename Int, typename... Ints>
auto ITensor::
cplx(Int iv1, Ints... ivs) const
    -> stdx::enable_if_t<std::is_integral<Int>::value && stdx::and_<std::is_integral<Ints>...>::value,Cplx>
    {
    if(!store()) Error("tensor storage unallocated");

    constexpr size_t size = sizeof...(ivs)+1;
    auto ints = std::array<Int,size>{{iv1,static_cast<int>(ivs)...}};
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
    //using index_type = typename IVal::index_type;
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
primeLevel(ITensor A, 
           VarArgs&&... vargs)
    {
    A.primeLevel(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
primeExcept(ITensor A, 
            VarArgs&&... vargs)
    {
    A.primeExcept(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
noprime(ITensor A, 
        VarArgs&&... vargs)
    {
    A.noprime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
mapprime(ITensor A, 
         VarArgs&&... vargs)
    {
    A.mapprime(std::forward<VarArgs>(vargs)...);
    return A;
    }

template<typename... VarArgs>
ITensor
sim(ITensor A, 
    VarArgs&&... vargs)
    {
    A.sim(std::forward<VarArgs>(vargs)...);
    return A;
    }


template<typename Cond>
Index
findindex(ITensor const& T, 
          Cond && cond)
    {
    for(auto& i : T.inds())
        if(cond(i)) return i;
    return Index{};
    }

//Find index of tensor A (of optional type t) 
//which is shared with tensor B
Index
commonIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t);


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
          IndexType type);

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
reindex(ITensor const& cT, 
        Index o1, Index n1, 
        Inds... inds) 
    {
    constexpr size_t size = 2+sizeof...(inds);
    auto ipairs = std::array<Index,size>{{o1,n1,static_cast<Index>(inds)...}};

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

IndexSet inline
moveToFront(IndexSet const& isf, IndexSet const& is)
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

    auto iso = IndexSet(r);

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

IndexSet inline
moveToBack(IndexSet const& isb, IndexSet const& is)
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

    auto iso = IndexSet(r);

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
template <typename... Indxs>
auto ITensor::
order(Index const& ind1, Indxs const&... inds)
    -> stdx::enable_if_t<not stdx::and_<std::is_same<Index, Indxs>...>::value,ITensor&>
    {
    static constexpr auto size = 1+(sizeof...(inds)-1);
    auto isf = IndexSet(size);
    detail::getDotInds(isf.begin(),ind1,std::forward<Indxs const&>(inds)...);
    order(detail::moveToFront(isf,this->inds()));
    return *this;
    }

//Version of order accepting syntax: T.order(i,j,k)
template <typename... Indxs>
auto ITensor::
order(Index const& ind1, Indxs const&... inds)
    -> stdx::enable_if_t<stdx::and_<std::is_same<Index, Indxs>...>::value,ITensor&>
    {
    order(IndexSet(ind1, inds...));
    return *this;
    }

//Version of order accepting syntax: T.order("...",j,k)
template <typename... Indxs>
ITensor& ITensor::
order(std::string const& dots, Indxs const&... inds)
    {
    if(dots != "...")
        Error(format("Wrong string passed to order (expected '...', got '%s')",dots));
    order(detail::moveToBack(IndexSet(inds...),this->inds()));
    return *this;
    }

//order function which returns a new ITensor 
template<typename... Indxs>
ITensor
order(ITensor A, Indxs const&... inds)
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
    auto len = minM(is);
    return ITensor(std::move(is),DiagReal(len,1.));
    }

template<typename Container, typename... Inds, class>
ITensor
diagTensor(Container const& C, 
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


} // namespace itensor


#endif
