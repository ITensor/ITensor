//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SKIP_ITERATOR_H
#define __ITENSOR_SKIP_ITERATOR_H

#include <iostream>
#include <iterator>

namespace itensor {
namespace detail {

// template<typename Iterator>
// struct ValidCheck
//     {
//     bool
//     operator()(const Iterator& it) const { return !it->valid(); }
//     };


template<typename Iterator,
         class Skip>
class SkipIterator
    {
    private:

    typedef std::iterator_traits<Iterator> __traits_type;

    public:

    typedef typename __traits_type::iterator_category iterator_category;
    typedef typename __traits_type::value_type value_type;
    typedef typename __traits_type::difference_type difference_type;
    typedef typename __traits_type::reference reference;
    typedef typename __traits_type::pointer pointer;

    typedef size_t size_type;

    public:

    //
    //  constructors
    //

    SkipIterator()
        { }

    ~SkipIterator()
        { }

    SkipIterator(Iterator start, Iterator end)
        : 
        curr_(start), 
        end_(end)
        {
        initialize_();
        }

    // convert from another SkipIterator with
    // compatible iterator type
    template<typename Iter>
    SkipIterator(const SkipIterator<Iter,Skip>& x)
        : 
        curr_(x.curr_), 
        end_(x.end_)
        { 
        }

    bool
    valid() const
        {
        return curr_ != end_;
        }

    //
    // comparison: general iterator requirements
    //

    bool 
    operator==(const SkipIterator& x) const
        {
        return curr_ == x.curr_;
        }

    bool 
    operator!=(const SkipIterator& x) const
        {
        return curr_ != x.curr_;
        }

    //
    // comparison: random access iterator requirements
    //

    bool 
    operator<(const SkipIterator& x) const
        {
        return curr_ < x.curr_;
        }

    bool operator<=(const SkipIterator& x) const
        {
        return curr_ <= x.curr_;
        }

    bool 
    operator>(const SkipIterator& x) const
        {
        return curr_ > x.curr_;
        }

    bool 
    operator>=(const SkipIterator& x) const
        {
        return curr_ >= x.curr_;
        }

    //
    //  access: forward iterator requirements
    //

    reference 
    operator*() const
        {
        return *curr_;
        }

    Iterator 
    operator->() const
        {
        return curr_;
        }

    SkipIterator& 
    operator++()
        {
        increment_();
        return *this;
        }

    SkipIterator  
    operator++(int)
        {
        SkipIterator save(*this);
        increment_();
        return save;
        }

    void
    swap(SkipIterator& x)
        {
        std::swap(curr_,x.curr_);
        std::swap(end_,x.end_);
        }

    private:

    void 
    increment_()
        {
        while(valid())
            {
            ++curr_;
            if(!valid() || !Skip()(curr_)) return;
            }
        }

    void 
    initialize_()
        {
        if(valid())
            {
            if(Skip()(curr_)) increment_();
            }
        }

    private:

    //
    //  member variables
    //

    /// iterator to current (keep to make fast access)
    Iterator curr_;

    /// iterator to current (keep to make fast access)
    Iterator end_;

    template <typename _i, typename _c>
    friend class SkipIterator;

    };

template<typename Iterator, class Skip>
std::ostream&
operator<<(std::ostream& s, const SkipIterator<Iterator,Skip>& si)
    {
    return s << "SkipIterator";
    }

} //namespace itensor::detail
} //namespace itensor

#endif 
