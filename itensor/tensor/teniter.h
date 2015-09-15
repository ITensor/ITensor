//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TENITER_H_
#define __ITENSOR_TENITER_H_

#include "itensor/matrix/rangeiter.h"

namespace itensor {

template<class Ptr, class RangeT>
class TenIter
    { 
    public:
    using value_type = typename std::iterator_traits<Ptr>::value_type;
    using reference = typename std::iterator_traits<Ptr>::reference;
    using difference_type = typename std::iterator_traits<Ptr>::difference_type;
    using pointer = typename std::iterator_traits<Ptr>::pointer;
    using iterator_category = std::forward_iterator_tag;
    using range_type = RangeT;
    using range_iter = typename range_type::const_iterator;
    using storage_type = DataRange<typename std::remove_pointer<Ptr>::type>;
    private:
    storage_type d_;
    range_iter it_;
    public: 

    TenIter() { }

    TenIter(storage_type d, range_type const& r) 
      : d_(d), 
        it_(r.begin()) 
        { }  

    reference 
    operator*() 
        { 
        return d_[it_.offset()];
        }

    TenIter& 
    operator++() { increment(); return *this; } 

    TenIter 
    operator++(int) { auto ct = *this; ct.increment(); return ct; } 

    //pointer
    //data() const { return p_; }

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
        end.it_ = r.end();
        return end;
        }
    }; 

} //namespace itensor

#endif

