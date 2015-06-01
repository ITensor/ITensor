//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATRIX_STRIDE_ITER_H_
#define __ITENSOR_MATRIX_STRIDE_ITER_H_

#include <iterator> 

namespace itensor {

//stride_iter based on the example by Cogswell, Diggins, and Stevens

template<class T> 
class stride_iter
    { 
    public:
    using value_type = typename std::iterator_traits<T>::value_type;
    using reference = typename std::iterator_traits<T>::reference;
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using iterator_category = std::random_access_iterator_tag;
    private:
    pointer p_; 
    long stride_; 
    public: 

    stride_iter() : p_(nullptr), stride_(0) { }

    stride_iter(const stride_iter& other) : p_(other.p_), stride_(other.stride_) { } 

    stride_iter(pointer p, long stride) : p_(p), stride_(stride) { }  

    pointer
    data() const { return p_; }
    long
    stride() const { return stride_; }

    stride_iter& 
    operator++() { p_ += stride_; return *this; } 
    stride_iter 
    operator++(int) { auto tmp = *this; p_ += stride_; return tmp; } 
    stride_iter& 
    operator+=(difference_type x) { p_ += x * stride_; return *this; } 
    stride_iter& 
    operator--( ) { p_ -= stride_; return *this; } 
    stride_iter 
    operator--(int) { auto tmp = *this; p_ -= stride_; return tmp; } 
    stride_iter& 
    operator-=(difference_type x) { p_ -= x * stride_; return *this; } 
    reference 
    operator[](difference_type n) { return p_[n * stride_]; } 
    reference 
    operator*() { return *p_; }  
    }; 

template <typename T>
bool 
operator==(const stride_iter<T>& x, const stride_iter<T>& y) 
    { assert(x.stride() == y.stride()); return x.data() == y.data(); } 
template <typename T>
bool 
operator!=(const stride_iter<T>& x, const stride_iter<T>& y) 
    { assert(x.stride() == y.stride()); return x.data() != y.data(); } 
template <typename T>
bool 
operator<(const stride_iter<T>& x, const stride_iter<T>& y) 
    { assert(x.stride() == y.stride()); return x.data() < y.data(); } 
template <typename T>
typename stride_iter<T>::difference_type 
operator-(const stride_iter<T>& x, const stride_iter<T>& y) 
    { assert(x.stride() == y.stride()); return (x.data() - y.data()) / x.stride(); } 
template <typename T>
stride_iter<T> 
operator+(const stride_iter<T>& x, typename stride_iter<T>::difference_type y) 
    { assert(x.stride() == y.stride()); return x += y * x.stride(); } 
template <typename T>
stride_iter<T> 
operator+(typename stride_iter<T>::difference_type x, const stride_iter<T>& y) 
    { assert(x.stride() == y.stride()); return y += x * x.stride(); } 

} //namespace itensor

#endif
