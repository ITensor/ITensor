//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_COUNT_H_
#define __ITENSOR_COUNT_H_

namespace itensor {

template <typename T>
struct CountHelper
    {
    constexpr
    CountHelper(const T& b,
                const T& e)
        : 
        curr_(b),
        end_(e)
        { }

    const T& 
    operator*() const { return curr_; }

    CountHelper& 
    operator++() 
        { 
        ++curr_; 
        return *this;
        }

    CountHelper 
    begin() const { return CountHelper(curr_,end_); }
    CountHelper 
    end() const { return CountHelper(end_,end_); }

    bool
    operator!=(const CountHelper& other)
        {
        return curr_ != other.curr_;
        }

    private:
    T curr_;
    T end_;
    };

template <typename T> constexpr
CountHelper<T>
count(T end)
    {
    return CountHelper<T>(0,end);
    }

template <typename ST, typename T> constexpr
CountHelper<T>
count(ST start, T end)
    {
    return CountHelper<T>(T(start),end);
    }
 
template <typename C> constexpr
auto
index(const C& container) -> CountHelper<decltype(container.size())>
    {
    using size_type = decltype(container.size());
    return CountHelper<size_type>(0,container.size());
    }

};

#endif
