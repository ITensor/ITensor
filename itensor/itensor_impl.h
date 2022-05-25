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
#ifndef __ITENSOR_ITENSOR_INTERFACE_IMPL_H_
#define __ITENSOR_ITENSOR_INTERFACE_IMPL_H_

#include "itensor/itdata/task_types.h"
#include "itensor/tensor/contract.h"

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

template <typename... Inds>
ITensor::
ITensor(Index  const& i1,
        Inds const&... inds)
  : is_(i1,inds...)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

template <typename... Inds>
ITensor::
ITensor(QN q, Index  const& i1,
        Inds const&... inds)
    {
    *this = ITensor(q,IndexSet(i1,inds...));
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
eltC(IV const& iv1, IVs&&... ivs) const
    -> stdx::if_compiles_return<Cplx,decltype(iv1.index),decltype(iv1.val)>
    {
    constexpr size_t size = sizeof...(ivs)+1;
    auto vals = std::array<IndexVal,size>{{static_cast<IndexVal>(iv1),
                                           static_cast<IndexVal>(ivs)...}};
    if(size != size_t(inds().order()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : vals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to elt/eltC (expected %d, got %d)",
                     inds().order(),size));
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
        println("too big for real in eltC(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in eltC(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }


template<typename Int>
auto ITensor::
eltC(std::vector<Int> const& ints) const
    -> stdx::enable_if_t<std::is_integral<Int>::value,Cplx>
    {
    if(!store()) Error("tensor storage unallocated");

    auto size = ints.size();
    if(size != size_t(inds().order()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        print("Indices provided = ");
        for(auto i : ints) print(" ",i);
        println("\n---------------------------------------------");
        Error(format("Wrong number of ints passed to elt/eltC (expected %d, got %d)",
                     inds().order(),size));
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
        println("too big for real in eltC(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in eltC(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

template<typename Int, typename... Ints>
auto ITensor::
eltC(Int iv1, Ints... ivs) const
    -> stdx::enable_if_t<std::is_integral<Int>::value 
                     && stdx::and_<std::is_integral<Ints>...>::value,Cplx>
    {
    return this->eltC(std::vector<Int>{{iv1,static_cast<int>(ivs)...}});
    }

template <typename... IVals>
Real ITensor::
elt(IVals&&... ivs) const
    {
    if(itensor::isComplex(*this)) Error("Cannot call .elt(...) on an ITensor with complex storage. Please use .eltC(...) instead");
    //TODO: make a specialized elt(...) version
    auto z = eltC(std::forward<IVals>(ivs)...);
    //if(fabs(z.imag()) > 1E-15 && fabs(z.imag()) > 1E-14*fabs(z.real()))
    //    {
    //    printfln("element = (%.5E,%.5E)",z.real(),z.imag());
    //    //Error("tensor is Complex valued, use .eltC(...) method");
    //    throw ITError("tensor is complex valued, use .eltC(...) method");
    //    }
    return z.real();
    }

template <typename... IVals>
Real ITensor::
real(IVals&&... ivs) const
    {
    return elt(std::forward<IVals>(ivs)...);
    }

template <typename... IVals>
Cplx ITensor::
cplx(IVals&&... ivs) const
    {
    return eltC(std::forward<IVals>(ivs)...);
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

template<typename Ints>
void
checkEltFluxInts(ITensor const& A, Ints const& ints)
    {
    if(hasQNs(A))
      {
      QN elt_flux;
      auto indsA = inds(A);
      for(auto i : range1(order(A)))
          {
          auto iv = indsA(i)(ints[i-1]+1);
          elt_flux += dir(iv)*qn(iv);
          }
      if(elt_flux != flux(A))
          {
          println("Trying to set element: ");
          for(auto i : range1(order(A)))
            println("Index: ", indsA(i), ", Val: ",ints[i-1]+1);
          println("Element flux is: ",elt_flux);
          println("ITensor flux is: ",flux(A));
          Error("In .set, cannot set element with flux different from ITensor flux");
          }
      }
    return;
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
    if(size != size_t(inds().order())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : vals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().order(),size));
        }
    auto inds = IntArray(is_.order(),0);
    detail::permute_map(is_,vals,inds,
                        [](IndexVal const& iv) { return iv.val-1; });
    //TODO: if !store_ and !is_real, call allocCplx instead
    //and move this line after check for is_real
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    detail::checkEltFluxInts(*this,inds);
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
    if(size != size_t(inds().order())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        print("Indices provided =");
        for(auto& i : ints) print(" ",1+i);
        println();
        println("---------------------------------------------");
        Error(format("Wrong number of ints passed to set (expected %d, got %d)",
                     inds().order(),size));
        }
    //TODO: if !store_ and !is_real, call allocCplx instead
    //and move this line after check for is_real
    if(!store_) detail::allocReal(*this,ints);
    scaleTo(1.);
    detail::checkEltFluxInts(*this,ints);
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
    fixBlockDeficient();
    scaleTo(1);
    doTask(GenerateIT<decltype(f)>{std::forward<Func>(f)},store_);
    return *this;
    }

template <typename Func>
ITensor& ITensor::
apply(Func&& f)
    {
    fixBlockDeficient();
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

namespace detail {

    template <typename... IVals>
    ITensor
    IndexValsToITensor(IndexVal const& iv1,
                       IVals const&... rest)
        {
        const constexpr auto size = 1+sizeof...(rest);
        auto ivs = stdx::make_array(iv1,rest...);
        //TODO: try directly making inds as iv1.index,(rest.index)...
        auto inds = std::array<Index,size>{};
        for(size_t j = 0; j < size; ++j) inds[j] = ivs[j].index;
        auto D = ITensor{IndexSet(inds)};
        return D;
        }

} //namespace detail

template <typename... IVals>
ITensor
setElt(Real el,
       IndexVal const& iv1, 
       IVals const&... rest)
    {
    auto D = detail::IndexValsToITensor(iv1, rest...);
    D.set(iv1,rest...,el);
    return D;
    }

template <typename... IVals>
ITensor
setElt(Cplx el,
       IndexVal const& iv1, 
       IVals const&... rest)
    {
    auto D = detail::IndexValsToITensor(iv1, rest...);
    D.set(iv1,rest...,el);
    return D;
    }

template <typename... IVals>
ITensor
setElt(IndexVal const& iv1, 
       IVals const&... rest)
    {
    return setElt(1.,iv1,rest...);
    }

template<typename... VarArgs>
Real
elt(ITensor A, 
    VarArgs&&... vargs)
    {
    return A.elt(std::forward<VarArgs>(vargs)...);
    }

template<typename... VarArgs>
Cplx
eltC(ITensor A, 
     VarArgs&&... vargs)
    {
    return A.eltC(std::forward<VarArgs>(vargs)...);
    }

//
//  I was unable to get these to instance with a varargs parameter pack.
//
template <> inline
Real
eltT(ITensor A)
    {
        return A.elt();
    }   

template <> inline
Complex
eltT(ITensor A)
    {
        return A.eltC();
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

template<typename T, typename... CtrArgs>
ITensor::storage_ptr
readType(std::istream& s, CtrArgs&&... args)
    {
    T t(std::forward<CtrArgs>(args)...);
    read(s,t);
    return newITData<T>(std::move(t));
    }

#ifdef ITENSOR_USE_HDF5

template<typename T>
ITensor::storage_ptr
h5_readStore(h5::group g, std::string const& name)
    {
    return newITData<T>(h5_read<T>(g,name));
    }

#endif

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

#ifdef ITENSOR_USE_HDF5

struct H5Write
    {
    h5::group& parent;
    std::string name;

    H5Write(h5::group& parent_,std::string const& name_) : parent(parent_), name(name_) { }
    };
inline const char*
typeNameOf(H5Write const&) { return "H5Write"; }

template<typename D>
auto
doTask(H5Write & W, D const& d)
    -> stdx::if_compiles_return<void,decltype(itensor::h5_write(W.parent,W.name,d))>
    {
    h5_write(W.parent,W.name,d);
    }

#endif 

template<typename Container, class>
ITensor
diagITensor(Container const& C, 
            IndexSet const& is)
    { 
    if( not hasQNs(is) )
      {
#ifdef DEBUG
      using size_type = decltype(C.size());
      //Compute min of all index dimensions
      auto mindim = dim(is[0]);
      for(const auto& ind : is)
          if(dim(ind) < mindim) mindim = dim(ind);
      if(C.size() != size_type(mindim))
          {
          println("mindim = ",mindim);
          println("C.size() = ",C.size());
          Error("Wrong size of data in diagonal ITensor constructor");
          }
#endif
      using value_type = typename Container::value_type;
      return ITensor(std::move(is),Diag<value_type>(C.begin(),C.end()));
      }
    else
      {
      Error("diagITensor constructor not yet implemented for QNs");
      return ITensor();
      }
    }

template<typename V>
TenRef<Range,V>
getBlock(ITensor & T,
         Block block_ind)
    {
    if(block_ind.size() != size_t(T.order())) Error("Mismatched number of indices and ITensor order");
    if(not T.store())
        {
        QN q;
        for(auto n : range(block_ind))
            {
            auto& I = inds(T)[n];
            q += qn(I,block_ind[n])*dir(I);
            }
        T = ITensor(inds(T),QDense<V>(inds(T),q));
        }
    //Interface is 1-indexed; switch to 0-indexed
    for(auto& i : block_ind) { i -= 1; }
    auto G = GetBlock<V>(inds(T),block_ind);
    return doTask(G,T.store());
    }

//
// Deprecated
//

template<typename Container, typename... Inds, class>
ITensor
diagTensor(Container const& C, 
           Index const& i1,
           Inds&&... inds)
    { 
    Global::warnDeprecated("diagTensor(Container,Index,...) is deprecated in favor of diagITensor(Container,Index,...)");
    return diagITensor(C,IndexSet(i1,std::forward<Inds>(inds)...));
    }

template <typename... Inds>
ITensor
randomTensor(Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensor(Index,...) is deprecated in favor of randomITensor(Index,...)");
    return randomITensor(i1,std::forward<Inds>(inds)...);
    }
template <typename... Inds>
ITensor
randomTensorC(Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensorC(Index,...) is deprecated in favor of randomITensorC(Index,...)");
    return randomITensorC(i1,std::forward<Inds>(inds)...);
    }

template <typename... Inds>
ITensor
randomTensor(QN q, Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensor(QN,Index,...) is deprecated in favor of randomITensor(QN,Index,...)");
    return randomITensor(q,i1,std::forward<Inds>(inds)...);
    }
template <typename... Inds>
ITensor
randomTensorC(QN q, Index const& i1, Inds&&... inds)
    {
    Global::warnDeprecated("randomTensorC(QN,Index,...) is deprecated in favor of randomITensorC(QN,Index,...)");
    return randomITensorC(q,i1,std::forward<Inds>(inds)...);
    }

template<typename... Inds>
ITensor
reindex(ITensor const& cT, 
        Index o1, Index n1, 
        Inds... indxs)
    {
    Error("Error: reindex(ITensor,Index,Index,...) is deprecated in favor of replaceInds(ITensor,Index,Index,...)");
    return ITensor();
    }

} // namespace itensor


#endif
