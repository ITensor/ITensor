//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_COUNT_H_
#define __ITENSOR_COUNT_H_

namespace itensor {

namespace detail {

template <typename size_type>
class CountHelper
    {
    size_type curr_;
    size_type end_;
    public:

    constexpr
    CountHelper(size_type b,
                size_type e)
        : 
        curr_(b),
        end_(e)
        { }

    const size_type&
    operator*() const { return curr_; }

    CountHelper& 
    operator++() 
        { 
        ++curr_; 
        return *this;
        }

    bool
    operator!=(const CountHelper& other) const
        {
        return curr_ != other.curr_;
        }

    CountHelper 
    begin() const { return CountHelper(curr_,end_); }

    CountHelper 
    end() const { return CountHelper(end_,end_); }
    };

} //namespace detail

template <typename T> constexpr
detail::CountHelper<T>
count(T end)
    {
    return detail::CountHelper<T>(0,end);
    }

template <typename ST, typename T> constexpr
detail::CountHelper<T>
count(ST start, T end)
    {
    return detail::CountHelper<T>(T(start),end);
    }
 
template <typename C> constexpr
auto
index(const C& container) -> detail::CountHelper<decltype(container.size())>
    {
    using size_type = decltype(container.size());
    return detail::CountHelper<size_type>(0,container.size());
    }

template <typename T> constexpr
detail::CountHelper<T>
count1(T end)
    {
    return detail::CountHelper<T>(1,1+end);
    }

template <typename ST, typename T> constexpr
detail::CountHelper<T>
count1(ST start, T end)
    {
    return detail::CountHelper<T>(start,1+end);
    }

}

#endif
