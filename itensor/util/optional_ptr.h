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
#ifndef __ITENSOR_OPTIONAL_PTR_H
#define __ITENSOR_OPTIONAL_PTR_H 1

#include <memory>
#include <type_traits>

namespace itensor {

/**

  optional_ptr<T> functions either as a raw unmanaged pointer
  or as a smart pointer of type managed_ptr depending on whether
  it is initialized through the set_external method or
  the set_managed method


*/
template <typename T, typename managed_ptr = std::unique_ptr<T>>
class optional_ptr
    {
    public:

    optional_ptr() : p_(nullptr) { }

    optional_ptr(const optional_ptr& other,
                 typename std::enable_if<std::is_copy_constructible<managed_ptr>::value>::type* = 0)
        :
        p_(other.p_),
        mp_(other.mp_)
        { }

    optional_ptr(optional_ptr&& other)
        :
        p_(other.p_),
        mp_(std::move(other.mp_))
        { }

    template <typename std::enable_if<std::is_copy_constructible<managed_ptr>::value>::type* = nullptr>
    optional_ptr&
    operator=(const optional_ptr& other)
        { 
        p_ = other.p_;
        mp_ = other.mp_;
        }

    T&
    operator*() const { return *p_; }

    T*
    operator->() const { return p_; }

    void
    set_managed(T* new_p)
        {
        mp_ = std::move(managed_ptr(new_p));
        p_ = mp_.get();
        }

    void
    set_external(T* ext_p)
        {
        p_ = ext_p;
        mp_.reset();
        }

    private:
    T* p_;
    managed_ptr mp_;
    };

} // namespace itensor

#endif
