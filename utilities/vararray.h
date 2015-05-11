//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_VARARRAY_H
#define __ITENSOR_VARARRAY_H

#include "print.h"

#ifdef DEBUG
#define BOUNDS_CHECK(X) bounds_check(X);
#else
#define BOUNDS_CHECK(X)
#endif

#ifdef DEBUG
#define SIZE_CHECK check_size();
#else
#define SIZE_CHECK
#endif

#ifdef DEBUG
#define EMPTY_CHECK check_empty();
#else
#define EMPTY_CHECK
#endif

namespace itensor {

template<typename T, size_t MaxSize>
class VarArray
    {
    public:
    using storage_type = std::array<T,MaxSize>;
    using value_type = typename storage_type::value_type;
    using size_type = typename storage_type::size_type;
    using difference_type = typename storage_type::difference_type;
    using reference = typename storage_type::reference;
    using const_reference = typename storage_type::const_reference;
    using pointer = typename storage_type::pointer;
    using const_pointer = typename storage_type::const_pointer;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;
    using reverse_iterator = typename storage_type::reverse_iterator;
    using const_reverse_iterator = typename storage_type::const_reverse_iterator;
    private:
    storage_type store_;
    size_t size_;
    public:

    VarArray() : size_(0) { }

    VarArray(size_t size) : size_(size) { }

    VarArray(size_t size,
              const_reference value) 
        : 
        size_(size) 
        { 
        store_.fill(value);
        }

    size_t
    size() const { return size_; }

    size_t constexpr
    max_size() const { return MaxSize; }

    void
    resize(size_t new_size) { size_ = new_size; }

    void
    clear() { size_ = 0; }

    void
    push_back(const_reference val) { store_[size_] = val; ++size_; SIZE_CHECK }

    void
    push_back(value_type&& val) { store_[size_] = std::move(val); ++size_; SIZE_CHECK }

    void
    assign(size_t count, const_reference val) { size_=count; store_.fill(val); SIZE_CHECK }

    explicit operator bool() const { return bool(size_); }

    reference
    operator[](size_t i) { BOUNDS_CHECK(i) return store_[i]; }

    const_reference
    operator[](size_t i) const { BOUNDS_CHECK(i) return store_[i]; }

    reference
    at(size_t i) { BOUNDS_CHECK(i) return store_[i]; }

    const_reference
    at(size_t i) const { BOUNDS_CHECK(i) return store_[i]; }

    reference
    front() { EMPTY_CHECK return store_.front(); }

    const_reference
    front() const { EMPTY_CHECK return store_.front(); }

    reference
    back() { EMPTY_CHECK return store_[size_==0 ? 0 : size_-1]; }

    const_reference
    back() const { EMPTY_CHECK return store_[size_==0 ? 0 : size_-1]; }

    pointer
    data() { EMPTY_CHECK return store_[0]; }

    const_pointer
    data() const { EMPTY_CHECK return store_[0]; }

    bool
    empty() const { return 0==size_; }

    void
    fill(const_reference val) { store_.fill(val); }

    void
    swap(VarArray& other) { std::swap(size_,other.size_); store_.swap(other.store_); }

    iterator
    begin() { return store_.data(); }

    iterator
    end() { return store_.data()+size_; }

    const_iterator
    begin() const { return store_.data(); }

    const_iterator
    end() const { return store_.data()+size_; }

    const_iterator
    cbegin() const { return store_.data(); }

    const_iterator
    cend() const { return store_.data()+size_; }

    private:
    void
    bounds_check(size_t i) const
        {
        if(i >= size_) throw std::out_of_range("index out of range in VarArray");
        }
    void
    check_size() const
        {
        if(size_==MaxSize) throw std::runtime_error("VarArray overflow, increase MaxSize");
        }
    void
    check_empty() const
        {
        if(size_==0) throw std::runtime_error("VarArray is empty");
        }
    };

#undef EMPTY_CHECK
#undef SIZE_CHECK
#undef BOUNDS_CHECK

} //namespace itensor

#endif
