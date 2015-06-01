//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INDEXSET_ITER_H
#define __ITENSOR_INDEXSET_ITER_H

#include <iterator> 

namespace itensor {

template<typename index_type_> 
class IndexSetIter
    { 
    public:
    using index_type = std::remove_const_t<index_type_>;
    using value_type = index_type;
    using reference = index_type_&;
    using difference_type = std::ptrdiff_t;
    using pointer = index_type_*;
    using iterator_category = std::random_access_iterator_tag;
    using range_ptr = typename RangeT<index_type>::value_type*;
    using const_range_ptr = const typename RangeT<index_type>::value_type*;
    using data_ptr = std::conditional_t<std::is_const<index_type_>::value,
                                      const_range_ptr,
                                      range_ptr>;
    private:
    data_ptr p_; 
    public: 

    IndexSetIter() : p_(nullptr) { }

    explicit
    IndexSetIter(data_ptr p) : p_(p) { }

    const data_ptr
    data() const { return p_; }

    IndexSetIter& 
    operator++() 
        { 
        ++p_; 
        return *this; 
        } 

    IndexSetIter 
    operator++(int) 
        { 
        auto tmp = *this; //save copy of this
        ++p_; 
        return tmp; 
        } 

    IndexSetIter& 
    operator+=(difference_type x) 
        { 
        p_ += x;
        return *this; 
        } 

    IndexSetIter& 
    operator--( ) 
        { 
        --p_;
        return *this; 
        } 

    IndexSetIter 
    operator--(int) 
        { 
        auto tmp = *this; //save copy of this
        --p_;
        return tmp; 
        } 

    IndexSetIter& 
    operator-=(difference_type x) 
        { 
        p_ -= x;
        return *this; 
        } 

    reference 
    operator[](difference_type n) { return p_[n].ext; } 

    reference 
    operator*() { return p_->ext; }  
    }; 

template <typename T>
bool 
operator==(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.data() == y.data(); 
    } 

template <typename T>
bool 
operator!=(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.data() != y.data(); 
    } 

template <typename T>
bool 
operator<(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.data() < y.data(); 
    } 

template <typename T>
typename IndexSetIter<T>::difference_type 
operator-(const IndexSetIter<T>& x, const IndexSetIter<T>& y) 
    { 
    return x.data() - y.data();
    } 

template <typename T>
IndexSetIter<T> 
operator+(const IndexSetIter<T>& x, typename IndexSetIter<T>::difference_type d) 
    { 
    return x += d;
    } 

template <typename T>
IndexSetIter<T> 
operator+(typename IndexSetIter<T>::difference_type d, const IndexSetIter<T>& x) 
    { 
    return x += d;
    } 

} //namespace itensor

#endif
