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
#ifndef __ITENSOR_INDEX_IMPL_H
#define __ITENSOR_INDEX_IMPL_H

namespace itensor {

namespace detail {

struct ArrowM
    {
    Arrow dir = Neither;
    TagSet tags;
    long m = 0l;
    ArrowM(Arrow d, 
           TagSet t,
           long m_) : dir(d), tags(t), m(m_) { }
    };

ArrowM inline
fill(std::vector<QNInt> const& v,
     Arrow dir,
     TagSet tags)
    { 
    return ArrowM(dir,tags,0l);
    }

ArrowM inline
fill(std::vector<QNInt> const& v,
     Arrow dir)
    { 
    return ArrowM(dir,TagSet(),0l);
    }

ArrowM inline
fill(std::vector<QNInt> const& v,
     TagSet tags)
    { 
    return ArrowM(Out,tags,0l);
    }

ArrowM inline
fill(std::vector<QNInt> const& v)
    { 
    return ArrowM(Out,TagSet(),0l);
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
Index(QN const& q1, long size1,
      QN_Sizes const&... qnsizes)
    { 
    constexpr auto size = 1+sizeof...(qnsizes)/2;
    auto qi = stdx::reserve_vector<std::pair<QN,long>>(size);
    auto am = detail::fill(qi,q1,size1,qnsizes...);
    auto I = Index(am.m,am.tags);
    operator=(I);
    setDir(am.dir);
    makeStorage(std::move(qi));
    }

} //namespace itensor


#endif
