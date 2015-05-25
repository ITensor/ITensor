//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_INFARRAY_H
#define __ITENSOR_INFARRAY_H

#include <array>
#include <vector>
#include <iterator> 

#ifdef DEBUG
#define CHECK_IND(X) check_ind(X);
#else
#define CHECK_IND(X)
#endif

#ifdef DEBUG
#define CHECK_EMPTY check_empty();
#else
#define CHECK_EMPTY
#endif

namespace itensor {

template<typename T, size_t ArrSize>
class infarray_iter;

template<typename T, size_t ArrSize>
class InfArray
    {
    public:
    using storage_type = std::array<T,ArrSize>;
    using value_type = typename storage_type::value_type;
    using size_type = typename storage_type::size_type;
    using difference_type = typename storage_type::difference_type;
    using reference = typename storage_type::reference;
    using const_reference = typename storage_type::const_reference;
    using pointer = typename storage_type::pointer;
    using const_pointer = typename storage_type::const_pointer;
    using iterator = infarray_iter<T,ArrSize>;
    using const_iterator = infarray_iter<const T,ArrSize>;
    private:
    storage_type store_;
    size_t size_;
    std::vector<T> vec_;
    public:

    InfArray() : size_(0) { }

    InfArray(size_t size)
        {
        resize(size);
        }

    InfArray(size_t size,
             const_reference value) 
        { 
        assign(size,value);
        }

    InfArray(std::initializer_list<T> init) 
        { 
        resize(init.size());
        auto it = init.begin();
        for(auto& el : store_)
            {
            el = *it;
            ++it;
            }
        if(!vec_.empty())
            {
            for(auto& el : vec_)
                {
                el = *it;
                ++it;
                }
            }
        }

    size_t
    size() const { return size_; }

    size_t constexpr
    array_size() const { return ArrSize; }

    size_t
    vec_size() const { return vec_.size(); }

    void
    resize(size_t new_size) 
        { 
        size_ = new_size; 
        if(size_ <= ArrSize)
            vec_.clear();
        else
            vec_.resize(size_-ArrSize);
        }

    void
    clear() 
        { 
        size_ = 0; 
        vec_.clear();
        }

    void
    push_back(const_reference val) 
        { 
        if(size_ < ArrSize) store_[size_] = val; 
        else                vec_.push_back(val);
        ++size_; 
        }

    void
    push_back(value_type&& val) 
        { 
        if(size_ < ArrSize) store_[size_] = std::move(val); 
        else                vec_.emplace_back(std::move(val));
        ++size_; 
        }

    void
    assign(size_t count, 
           const_reference val) 
        { 
        resize(count);
        fill(val);
        }

    explicit operator bool() const { return bool(size_); }

    reference
    operator[](size_t i) { CHECK_IND(i) return i < ArrSize ? store_[i] : vec_[i-ArrSize]; }

    const_reference
    operator[](size_t i) const { CHECK_IND(i) return i < ArrSize ? store_[i] : vec_[i-ArrSize]; }

    reference
    at(size_t i) { return i < ArrSize ? store_[i] : vec_.at(i-ArrSize); }

    const_reference
    at(size_t i) const { return i < ArrSize ? store_[i] : vec_.at(i-ArrSize); }

    reference
    front() { CHECK_EMPTY return store_.front(); }

    const_reference
    front() const { CHECK_EMPTY return store_.front(); }

    reference
    back() { CHECK_EMPTY return vec_.empty() ? store_[size_-1] : vec_.back(); }

    const_reference
    back() const { CHECK_EMPTY return vec_.empty() ? store_[size_-1] : vec_.back(); }

    bool
    empty() const { return size_==0; }

    void
    fill(const_reference val) 
        { 
        store_.fill(val); 
        std::fill(vec_.begin(),vec_.end(),val);
        }

    void
    swap(InfArray& other) 
        { 
        std::swap(size_,other.size_); 
        store_.swap(other.store_); 
        vec_.swap(other.vec_);
        }

    iterator
    begin() 
        { 
        if(size_==0) 
            {
            return end();
            }
        else if(size_ < ArrSize)
            {
            return iterator(store_.data(),
                            store_.data()+size_,
                            vec_.data()); 
            }
        else
            {
            return iterator(store_.data(),
                            store_.data()+ArrSize,
                            vec_.data()); 
            }
        }

    iterator
    end() 
        { 
        return iterator(vec_.data()+vec_.size(),
                        store_.data()+size_,
                        vec_.data());
        }

    const_iterator
    begin() const
        { 
        if(size_==0) 
            {
            return end();
            }
        else if(size_ < ArrSize)
            {
            return const_iterator(store_.data(),
                                  store_.data()+size_,
                                  vec_.data()); 
            }
        else
            {
            return const_iterator(store_.data(),
                                  store_.data()+ArrSize,
                                  vec_.data()); 
            }
        }

    const_iterator
    end() const
        { 
        return const_iterator(vec_.data()+vec_.size(),
                              store_.data()+size_,
                              vec_.data());
        }

    const_iterator
    cbegin() const { return begin(); }

    const_iterator
    cend() const { return end(); }

    private:
    void
    check_ind(size_t i) const
        {
        if(i >= size_) throw std::out_of_range("index out of range in InfArray");
        }
    void
    check_empty() const
        {
        if(size_==0) throw std::runtime_error("InfArray is empty");
        }
    };

template<typename T, size_t ArrSize>
class infarray_iter
    {
    using parent = typename std::iterator<std::forward_iterator_tag, T>;
    public:
    using value_type = typename parent::value_type;
    using reference = typename parent::reference;
    using pointer = T*;
    using difference_type = typename parent::difference_type;
    using iterator_category = typename parent::iterator_category;
    private:
    pointer p_; 
    pointer array_end_;
    pointer vec_start_;
    public: 

    infarray_iter() : 
        p_(nullptr), 
        array_end_(nullptr), 
        vec_start_(nullptr) 
        { }

    infarray_iter(pointer start,
                  pointer array_end,
                  pointer vec_start) :
        p_(start),
        array_end_(array_end),
        vec_start_(vec_start)
        { }

    infarray_iter(const infarray_iter& other) : 
        p_(other.p_), 
        array_end_(other.array_end_), 
        vec_start_(other.vec_start_)
        { } 

    pointer
    data() const { return p_; }

    infarray_iter& 
    operator++() 
        { 
        ++p_;
        if(p_==array_end_) p_ = vec_start_;
        return *this; 
        } 
    infarray_iter 
    operator++(int) 
        { 
        auto tmp = *this; 
        operator++();
        return tmp; 
        } 
    reference 
    operator*() { return *p_; }  
    };

template <typename T, size_t ArrSize>
bool 
operator!=(const infarray_iter<T,ArrSize>& x, const infarray_iter<T,ArrSize>& y) 
    { return x.data() != y.data(); } 

} //namespace itensor

#undef CHECK_EMPTY
#undef CHECK_IND

#endif
