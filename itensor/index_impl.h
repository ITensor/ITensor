//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEX_IMPL_H
#define __ITENSOR_INDEX_IMPL_H

namespace itensor {

namespace detail {

struct ArrowM
    {
    Arrow dir = Neither;
    long m = 0l;
    ArrowM(Arrow d, long m_) : dir(d), m(m_) { }
    };

ArrowM inline
fill(std::vector<QNInt> const& v,
     Arrow dir = Out) 
    { 
    return ArrowM(dir,0l);
    }

template<typename... Rest>
ArrowM
fill(std::vector<QNInt> & v,
     QN const& q, 
     long size, 
     Rest const&... rest)
    {
    v.emplace_back(q,size);
    auto am = fill(v,rest...);
    am.m += size;
    return am;
    }

} //namespace detail

template<typename... QN_Sizes>
Index::
Index(std::string const& name, 
      QN const& q1, long size1,
      QN_Sizes const&... qnsizes)
    : type_(Link)
    { 
    constexpr auto size = 1+sizeof...(qnsizes)/2;
    auto qi = stdx::reserve_vector<std::pair<QN,long>>(size);
    auto am = detail::fill(qi,q1,size1,qnsizes...);
    dir_ = am.dir;
    auto I = Index(name,am.m);
    operator=(I);
    makeStorage(std::move(qi));
    }

} //namespace itensor


#endif
