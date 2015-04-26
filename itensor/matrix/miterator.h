//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MITERATOR_H_
#define __ITENSOR_MITERATOR_H_

#include <iterator> 
#include "mrange.h"

namespace itensor {

template<class T> 
class miterator
    { 
    public:
    using value_type = typename std::iterator_traits<T>::value_type;
    using reference = typename std::iterator_traits<T>::reference;
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using iterator_category = std::forward_iterator_tag;
    private:
    pointer p_; 
    long count_;
    mrange ind_; 
    public: 

    miterator() : p_(nullptr), count_(0) { }; 
    miterator(const miterator& other) : p_(other.p_), count_(other.count_), ind_(other.ind_) { } 
    miterator(pointer p, const mrange& ind) : p_(p), count_(0), ind_(ind) { }  

    pointer
    data() const { return p_; }
    const mrange&
    ind() const { return ind_; }

    miterator& 
    operator++() { increment(); return *this; } 
    miterator 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 
    reference 
    operator*() { return *p_; }  

    bool
    operator!=(const miterator& other) const { return count_!=other.count_; }
    bool
    operator==(const miterator& other) const { return count_==other.count_; }

    private:

    void
    increment()
        {
        p_ += ind_.rs;
        ++count_;
        if((count_%ind_.rn) == 0)
            {
            std::advance(p_,ind_.cs-ind_.rn*ind_.rs);
            }
        }
    public:
    //For developer use only; for making end iterator
    miterator(const mrange& ind) : p_(nullptr), count_(ind.area()), ind_(ind) { }
    }; 


template <typename T>
bool 
operator==(const miterator<T>& x, const miterator<T>& y) { assert(x.ind() == y.ind()); return x.data() == y.data(); } 
template <typename T>
bool 
operator!=(const miterator<T>& x, const miterator<T>& y) { assert(x.ind() == y.ind()); return x.data() != y.data(); } 

};

#endif
