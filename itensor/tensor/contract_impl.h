//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CONTRACT_IMPL_H__
#define __ITENSOR_CONTRACT_IMPL_H__

namespace itensor {

namespace detail {

//Helper class for ncprod.
//Holds an container of pointers to value_type,
//PtrInter interface pretends this is an actual 
//container of values (not pointers) instead
template<typename value_type_, size_t ArrSize>
class PtrInd
    {
    public:
    using value_type = value_type_;
    using pointer_type = value_type const*;
    using storage_type = InfArray<pointer_type,ArrSize>;
    using size_type = typename storage_type::size_type;
    private:
    storage_type ptrs_;
    public:

    PtrInd(size_type size)
      : ptrs_(size,nullptr)
        { }

    size_type
    size() const { return ptrs_.size(); }

    void
    set(size_type n, pointer_type p)
        {
        ptrs_[n] = p;
        }

    value_type const&
    operator[](size_type n) const
        {
        return *(ptrs_[n]);
        }
    };

} //namespace detail

//Non-contracting product
template<typename R, typename VA, typename VB>
void 
ncprod_impl(TenRefc<R,VA> A, Labels const& ai, 
            TenRefc<R,VB> B, Labels const& bi, 
            TenRef<R,common_type<VA,VB>>  C, Labels const& ci)
    {
    auto rA = rank(A),
         rB = rank(B),
         rC = rank(C);

    auto cb = rangeBegin(C.range());
    auto ce = rangeEnd(C.range());

    using value_type = stdx::remove_reference_t<decltype(cb[0])>;
    using PtrIndType = detail::PtrInd<value_type,Labels::arr_size()>;
    auto aind = PtrIndType(rA);
    auto bind = PtrIndType(rB);

    for(auto nc : range(rC))
        {
        for(auto na : range(rA))
            {
            if(ci[nc] == ai[na])
                {
                aind.set(na,&cb[nc]);
                break;
                }
            }
        for(auto nb : range(rB))
            {
            if(ci[nc] == bi[nb])
                {
                bind.set(nb,&cb[nc]);
                break;
                }
            }
        }

    auto pa = MAKE_SAFE_PTR(A.data(),A.size());
    auto pb = MAKE_SAFE_PTR(B.data(),B.size());
    auto pc = MAKE_SAFE_PTR(C.data(),C.size());
    for(; cb != ce; ++cb)
        {
        pc[cb.offset()] = pa[offset(A,aind)] * pb[offset(B,bind)];
        }
    }

template<class TA, class TB, class TC>
void 
ncprod(TA && A, Labels const& ai, 
       TB && B, Labels const& bi, 
       TC && C, Labels const& ci)
    {
    ncprod_impl(makeRef(A),ai,makeRef(B),bi,makeRef(C),ci);
    }

} // namespace itensor

#endif
