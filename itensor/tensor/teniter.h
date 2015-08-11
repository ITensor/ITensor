//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TENITER_H_
#define __ITENSOR_TENITER_H_

#include "rangeiter.h"

namespace itensor {

template<class T, class RangeT>
class TenIter
    { 
    public:
    using value_type = typename std::iterator_traits<T>::value_type;
    using reference = typename std::iterator_traits<T>::reference;
    using difference_type = typename std::iterator_traits<T>::difference_type;
    using pointer = typename std::iterator_traits<T>::pointer;
    using iterator_category = std::forward_iterator_tag;
    using range_type = RangeT;
    using range_iter = RangeIter<RangeT>;
    private:
    pointer p_; 
    range_iter it_;
    public: 

    TenIter() : p_(nullptr) { }

    TenIter(pointer p, range_type const& r) 
      : p_(p), 
        it_(r) 
        { }  

    reference 
    operator*() { return *(p_+it_.offset()); }  

    TenIter& 
    operator++() { increment(); return *this; } 

    TenIter 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 

    pointer
    data() const { return p_; }

    range_type const&
    range() const { return it_.range(); }

    bool
    notDone() const { return it_.notDone(); }

    bool
    operator!=(TenIter const& other) const { return it_!=other.it_; }

    bool
    operator==(TenIter const& other) const { return it_==other.it_; }

    private:

    void
    increment() { ++it_; }

    public:
    //For developer use only; for making end iterator
    TenIter static
    makeEnd(range_type const& r) 
        {
        TenIter end;
        end.it_ = range_iter::makeEnd(r);
        return end;
        }
    }; 

} //namespace itensor

#endif

