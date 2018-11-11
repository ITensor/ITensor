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
    IndexType type = Link;
    long m = 0l;
    ArrowM(Arrow d, 
           IndexType t,
           long m_) : dir(d), type(t), m(m_) { }
    };

ArrowM inline
fill(std::vector<QNInt> const& v,
     IndexType type = Link,
     Arrow dir = Out)
    { 
    return ArrowM(dir,type,0l);
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
    : type_(Link),
      dir_(Out)
    { 
    constexpr auto size = 1+sizeof...(qnsizes)/2;
    auto qi = stdx::reserve_vector<std::pair<QN,long>>(size);
    auto am = detail::fill(qi,q1,size1,qnsizes...);
    auto I = Index(name,am.m,am.type);
    operator=(I);
    dir(am.dir);
    makeStorage(std::move(qi));
    }

} //namespace itensor


#endif
