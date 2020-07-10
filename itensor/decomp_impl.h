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
#ifndef __ITENSOR_DECOMP_IMPL_H
#define __ITENSOR_DECOMP_IMPL_H

namespace itensor {

namespace detail {

Args inline
makeIndexSetArgs(std::vector<Index>& inds, Args const& args)
  {
  return args;
  }

Args inline
makeIndexSetArgs(std::vector<Index>& inds, Index const& i)
  {
  inds.push_back(i);
  return Args::global();
  }

template <typename... IndsArgs>
Args
makeIndexSetArgs(std::vector<Index>& inds,
                 Index const& i1,
                 IndsArgs&&... indsargs)
  {
  inds.push_back(i1);
  return makeIndexSetArgs(inds, indsargs...);
  }

template <typename... IndsArgs>
std::tuple<IndexSet,Args>
makeIndexSetArgs(Index const& i1, IndsArgs&&... indsargs)
  {
  auto inds = std::vector<Index>();
  auto args = makeIndexSetArgs(inds,i1,indsargs...);
  auto is = IndexSet(inds);
  return std::make_tuple(is,args);
  }

} //namespace detail

template <typename... IndsArgs>
std::tuple<ITensor,ITensor,ITensor>
svd(ITensor const& T, Index const& i1, IndsArgs&&... indsargs)
  {
  auto [is,args] = detail::makeIndexSetArgs(i1,indsargs...);
  return svd(T,is,args);
  }

template <typename... IndsArgs>
std::tuple<ITensor,ITensor>
qr(ITensor const& T, Index const& i1, IndsArgs&&... indsargs)
  {
  auto [is,args] = detail::makeIndexSetArgs(i1,indsargs...);
  return qr(T,is,args);
  }

template <typename... IndsArgs>
std::tuple<ITensor,ITensor>
factor(ITensor const& T, Index const& i1, IndsArgs&&... indsargs)
  {
  auto [is,args] = detail::makeIndexSetArgs(i1,indsargs...);
  return factor(T,is,args);
  }

} //namespace itensor

#endif
