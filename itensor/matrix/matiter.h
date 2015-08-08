//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MATITER_H_
#define __ITENSOR_MATITER_H_

#include <iterator> 
#include "mrange.h"

namespace itensor {

template<class T> 
class MatIter
    { 
    public:
    using value_type = typename std::iterator_traits<T>::value_type;
    using reference = typename std::iterator_traits<T>::reference;
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using iterator_category = std::forward_iterator_tag;
    using range_type = MRange;
    private:
    pointer p_; 
    long count_;
    range_type range_; 
    public: 

    MatIter() : p_(nullptr), count_(0) { }; 
    MatIter(MatIter const& other) : p_(other.p_), count_(other.count_), range_(other.range_) { } 
    MatIter(pointer p, range_type const& r) : p_(p), count_(0), range_(r) { }  

    pointer
    data() const { return p_; }

    range_type const&
    range() const { return range_; }

    MatIter& 
    operator++() { increment(); return *this; } 

    MatIter 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 

    reference 
    operator*() { return *p_; }  

    bool
    operator!=(MatIter const& other) const { return count_!=other.count_; }

    bool
    operator==(MatIter const& other) const { return count_==other.count_; }

    private:

    void
    increment()
        {
        p_ += range_.rs;
        ++count_;
        if((count_%range_.rn) == 0)
            std::advance(p_,range_.cs-range_.rn*range_.rs);
        }
    public:
    //For developer use only; for making end iterator
    MatIter(range_type const& ind) : p_(nullptr), count_(ind.area()), range_(ind) { }
    }; 


template <typename T>
bool 
operator==(MatIter<T> const& x, MatIter<T> const& y) { assert(x.range() == y.range()); return x.data() == y.data(); } 
template <typename T>
bool 
operator!=(MatIter<T> const& x, MatIter<T> const& y) { assert(x.range() == y.range()); return x.data() != y.data(); } 

}

#endif
