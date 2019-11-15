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
#ifndef __ITENSOR_VECTOR_NO_INIT_H
#define __ITENSOR_VECTOR_NO_INIT_H

#include <vector>

namespace itensor {

template <class T>
class uninitialized_allocator
  {
  public:
  typedef T value_type;

  uninitialized_allocator() noexcept { }

  template <class U>
  uninitialized_allocator(uninitialized_allocator<U> const&) noexcept { }

  T*
  allocate(std::size_t n)
    {
    return static_cast<T*>(::operator new(n * sizeof(T)));
    }

  void
  deallocate(T* p, std::size_t) noexcept
    {
    ::operator delete(static_cast<void*>(p));
    }

  template <class U>
  void construct(U*) noexcept
    {
    //TODO: do we want this trait check? It fails for std::complex<double>
    //static_assert(std::is_trivially_default_constructible<U>::value,
    //"This allocator can only be used with trivally default constructible types");
    }

  template <class U, class A0, class... Args>
  void
  construct(U* up, A0&& a0, Args&&... args) noexcept
    {
    ::new(up) U(std::forward<A0>(a0), std::forward<Args>(args)...);
    }

  bool
  operator==(uninitialized_allocator<T> const&) { return true; }

  bool
  operator!=(uninitialized_allocator<T> const&) { return false; }

  };

template<typename T>
using vector_no_init = std::vector<T,uninitialized_allocator<T>>;

} //namespace itensor

#endif
