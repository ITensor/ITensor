//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TEN_IMPL_H_
#define __ITENSOR_TEN_IMPL_H_

namespace itensor {

template<typename R,typename T>
TenRefc<R,T>& TenRefc<R,T>::
operator=(TenRefc const& t)
    {
    d_ = t.d_;
    if(t.ownRange())
        {
        range_ = t.range_;
        prange_ = &range_;
        }
    else
        {
        prange_ = t.prange_;
        }
    return *this;
    }

template<typename R, typename T>
TenRefc<R,T>& TenRefc<R,T>::
operator=(TenRefc && t)
    {
    d_ = t.d_;
    if(t.ownRange())
        {
        range_ = std::move(t.range_);
        prange_ = &range_;
        }
    else
        {
        prange_ = t.prange_;
        }
    return *this;
    }

template<typename R, typename T>
auto TenRefc<R,T>::
operator()() const -> reference
    { 
#ifdef DEBUG
    if(r() != 0) throw std::runtime_error("No indices passed to rank > 0 TenRef");
#endif
    return d_[0];
    }


template<typename R,typename T>
template <typename... Inds>
auto TenRefc<R,T>::
operator()(Inds&&... ii) const -> reference
    { 
    return d_[offset(*prange_,std::forward<Inds>(ii)...)]; 
    }

template<typename R,typename T>
auto TenRefc<R,T>::
operator[](size_t n) const -> reference
    { 
    return d_[n];
    }

template<typename R,typename T>
void TenRefc<R,T>::
pointTo(tensor_type const& t)
    {
    d_ = storage_type(t.data(),t.size());
    prange_ = &t.range();
    range_ = range_type{};
    }

template<typename R1, typename T1, 
         typename R2, typename T2>
void
checkCompatible(TenRefc<R1,T1> const& A, 
                TenRefc<R2,T2> const& B,
                std::string methodName = "")
    {
    auto methodstr = (methodName != "" ? format("in %s",methodName) : "");
    if(A.r() != B.r()) Error(format("Mismatched tensor ranks %s",methodstr));
    for(decltype(A.r()) n = 0; n < A.r(); ++n)
        if(A.extent(n) != B.extent(n))
            {
            printfln("A.extent(%d)=%d  B.extent(%d)=%d",n,A.extent(n),n,B.extent(n));
            Error(format("Mismatched tensor extent %s",methodstr));
            }
    }

template<typename R1, typename T1, 
         typename R2, typename T2, 
         typename Op>
void
transform(TenRefc<R1,T1> const& from, 
          TenRef<R2,T2>  const& to,
          Op&& op)
    {
#ifdef DEBUG
    checkCompatible(to,from,"transform");
#endif 
    using size_type = decltype(from.extent(0));
    auto r = to.r();
    if(r == 0)
        {
        op(*(from.data()),*(to.data()));
        return;
        }

    //find size and location of largest index of from
    size_type bigind = 0, 
              bigsize = from.extent(0);
    for(decltype(r) j = 1; j < r; ++j)
        if(bigsize < from.extent(j))
            {
            bigsize = from.extent(j); 
            bigind = j;
            }

    auto stepfrom = from.stride(bigind);
    auto stepto = to.stride(bigind);

    auto RB = RangeBuilder(r);
    for(decltype(r) i = 0; i < r; ++i)
        RB.setIndex(i,from.extent(i));
    //Leave bigind fixed to zero, will
    //increment manually in the loop below
    RB.setIndex(bigind,1);

    for(auto& i : RB.build())
        {
        //println("i = ",i);
        //printfln("to (offset,size) = (%d,%d)",offset(to,i),to.store().size());
        //printfln("from (offset,size) = (%d,%d)",offset(from,i),from.store().size());
        auto pto = MAKE_SAFE_PTR_OFFSET(to.data(),offset(to,i),to.store().size());
        auto pfrom = MAKE_SAFE_PTR_OFFSET(from.data(),offset(from,i),from.store().size());
        for(decltype(bigsize) b = 0; b < bigsize; ++b)
            {
            op(*pfrom,*pto);
            pto += stepto;
            pfrom += stepfrom;
            }
        }
    }

//Assign to referenced data
template<typename R1, typename R2, typename T>
void 
operator&=(TenRef<R1,T> const& A, TenRefc<R2,T> const& B)
    {
    transform(B,A,[](T b, T& a){ a = b; });
    }

//Assign to referenced data
template<typename R1, typename R2, typename T>
void
operator&=(TenRef<R1,T> const& A, Ten<R2,T> const& B)
    {
    transform(makeRef(B),A,[](T b, T& a){ a = b; });
    }

template<typename R1, typename R2,typename T>
void 
operator+=(TenRef<R1,T> const& A, TenRefc<R2,T> const& B)
    {
    transform(B,A,[](T b, T& a){ a += b; });
    }

template<typename R, typename T>
void
operator+=(TenRef<R,T> const& A, Ten<Range,T> const& B)
    {
    transform(makeRef(B),A,[](T b, T& a){ a += b; });
    }

template<typename R,typename T>
void
operator+=(Ten<Range,T> & A, TenRefc<R,T> const& B)
    {
    transform(B,makeRef(A),[](T b, T& a){ a += b; });
    }

//template<typename R, typename T>
//void
//operator+=(Ten<R,T> & A, Ten<R,T> const& B)
//    {
//    transform(makeRef(B),makeRef(A),[](T b, T& a){ a += b; });
//    }

template<typename R,typename T>
auto inline Ten<R,T>::
operator()() const -> value_type
    { 
#ifdef DEBUG
    if(r() != 0) throw std::runtime_error("No indices passed to rank > 0 Ten");
    if(data_.empty()) throw std::runtime_error("Empty storage in tensor when calling operator()");
#endif
    return data_.front();
    }

template<typename R,typename T>
template <typename... Inds>
auto Ten<R,T>::
operator()(Inds&&... ii) const -> const_reference
    { 
    return store()[offset(range_,std::forward<Inds>(ii)...)]; 
    }

template<typename R,typename T>
auto inline Ten<R,T>::
operator()() -> reference
    { 
#ifdef DEBUG
    if(r() != 0) throw std::runtime_error("No indices passed to rank > 0 Ten");
    if(data_.empty()) throw std::runtime_error("Empty storage in tensor when calling operator()");
#endif
    return data_.front(); 
    }

template<typename R,typename T>
auto inline Ten<R,T>::
operator[](size_t n) const -> const_reference
    { 
    return data_[n];
    }

template<typename R,typename T>
auto inline Ten<R,T>::
operator[](size_t n) -> reference
    { 
    return data_[n];
    }

template<typename R,typename T>
template <typename... Inds>
auto Ten<R,T>::
operator()(Inds&&... ii) -> reference
    { 
    return store()[offset(range_,std::forward<Inds>(ii)...)]; 
    }

template<typename R,typename T>
auto TenRef<R,T>::
operator[](size_t n) const -> reference
    { 
    return store()[n];
    }

template<typename R, typename V>
Real
norm(TenRefc<R,V> const& t)
    {
    if(isContiguous(t))
        {
        auto d = realData(t);
        return dnrm2_wrapper(d.size(),d.data());
        }
    Real nrm = 0;
    for(auto& el : t) nrm += std::norm(el);
    return std::sqrt(nrm);
    }

Tensor inline
scalarTen(Real val)
    {
    return Tensor{Tensor::storage_type(1,val),Range{}};
    }

namespace detail {

template<typename X>
auto
random() -> stdx::enable_if_t<std::is_same<X,Real>::value,Real>
    {
    return detail::quickran();
    }

template<typename X>
auto
random() -> stdx::enable_if_t<std::is_same<X,Cplx>::value,Cplx>
    {
    return Cplx(detail::quickran(),detail::quickran());
    }

}//namespace detail

template<typename R,typename V>
void
randomize(TenRef<R,V> const& t)
    {
    for(auto& el : t) el = detail::random<V>();
    }

template<typename R,typename V>
void
randomize(Ten<R,V> & t)
    {
    for(auto& el : t) el = detail::random<V>();
    }

template<typename R, typename V>
void
conjugate(TenRef<R,V> const& T)
    {
    if(isCplx<V>())
        {
        if(isContiguous(T))
            {
            auto de = T.data()+T.size();
            for(auto d = T.data(); d != de; ++d)
                {
                applyConj(*d);
                }
            }
        else
            {
            for(auto& el : T)
                {
                applyConj(el);
                }
            }
        }
    }

template<typename R, typename V>
void
conjugate(Ten<R,V> & T)
    {
    conjugate(makeRef(T));
    }

//return conjugated copy
template<typename R, typename V>
Ten<R,V>
conj(TenRefc<R,V> const& T)
    {
    auto CT = Ten<R,V>{T};
    conjugate(CT);
    return CT;
    }

template<typename R, typename V>
Ten<R,V>
conj(Ten<R,V> T)
    {
    conjugate(T);
    return T;
    }

template<typename R, typename V>
std::ostream&
printTensor(std::ostream & s, 
            TenRefc<R,V> const& T, 
            const char* typestr)
    {
    if(not T) return s << "(empty " << typestr << ")";
    if(T.r() == 0)
        {
        s << typestr << "\n() " << T();
        return s;
        }
    s << typestr << "\n";
    auto e = rangeEnd(T.range());
    for(auto i = rangeBegin(T.range()); i != e; ++i)
        {
        s << "(" << i.index() << ") " << formatVal(T(i)) << "\n";
        }
    return s;
    }

template<typename R,typename V>
std::ostream&
operator<<(std::ostream & s, TenRef<R,V> const& T)
    {
    return printTensor(s,T,"TenRef");
    }

template<typename R, typename V>
std::ostream&
operator<<(std::ostream & s, TenRefc<R,V> const& T)
    {
    return printTensor(s,T,"TenRefc");
    }

template<typename R, typename V>
std::ostream&
operator<<(std::ostream & s, Ten<R,V> const& T) 
    { 
    return printTensor(s,makeRefc(T),"Tensor");
    }

template<typename TenA, typename TenB,
         class = stdx::require<isTensor<TenA>,
                               isTensor<TenB>,
                               std::is_same<range_type<TenA>,range_type<TenB>>> >
auto
operator+(TenA && A, TenB && B)
    -> Ten<range_type<TenA>,common_type<TenA,TenB>>
    {
    using ResType = Ten<range_type<TenA>,common_type<TenA,TenB>>;
    ResType res;
    if(stdx::isRvalue<TenB>() || (isCplx(B) && isReal(A)))
        {
        res = ResType(std::forward<TenB>(B));
        operator+=(res,std::forward<TenA>(A));
        }
    else
        {
        res = ResType(std::forward<TenA>(A));
        operator+=(res,std::forward<TenB>(B));
        }
    return res;
    }

template<typename TenA, typename TenB,
         class = stdx::require<isTensor<TenA>,
                               isTensor<TenB>,
                               std::is_same<range_type<TenA>,range_type<TenB>>> >
auto
operator-(TenA && A, TenB && B)
    -> Ten<range_type<TenA>,common_type<TenA,TenB>>
    {
    using ResType = Ten<range_type<TenA>,common_type<TenA,TenB>>;
    ResType res;
    if(stdx::isRvalue<TenB>() || (isCplx(B) && isReal(A)))
        {
        res = ResType(std::forward<TenB>(B));
        res *= -1;
        operator+=(res,std::forward<TenA>(A));
        }
    else
        {
        res = ResType(std::forward<TenA>(A));
        operator-=(res,std::forward<TenB>(B));
        }
    return res;
    }

} //namespace itensor

#endif
