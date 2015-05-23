//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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
