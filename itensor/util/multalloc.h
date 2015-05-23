//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_MULTALLOC_H
#define __ITENSOR_MULTALLOC_H

#include <memory>

#ifdef DEBUG
#define CHECK_IND(X) check_ind(X);
#else
#define CHECK_IND(X)
#endif

#ifdef DEBUG
#define CHECK_SIZE(X) check_size(X);
#else
#define CHECK_SIZE(X)
#endif

#ifdef DEBUG
#define CHECK_NOT_FULL check_not_full();
#else
#define CHECK_NOT_FULL
#endif

#ifdef DEBUG
#define CHECK_NOT_EMPTY check_not_empty();
#else
#define CHECK_NOT_EMPTY
#endif

#ifdef DEBUG
#define CHECK_ALLOCATED check_allocated();
#else
#define CHECK_ALLOCATED
#endif

namespace itensor {

//
// //Sample usage:
// MultAlloc<Real,3> ma;
// ma.add(size0);
// ma.add(size1);
// ma.allocate(); //do single allocation for both memory ranges
// //ma.add(size2); //error: can't request more memory, idea is to only allocate once
// auto* p0 = ma[0];
// auto* p1 = ma[1];
//

template<typename T, size_t MaxNAlloc, typename Alloc = std::allocator<T>>
class MultAlloc
    {
    public:
    using size_type = std::size_t;
    using pointer = T*;
    private:
    struct SizeOff
        {
        size_type size = 0;
        size_type offset = 0;
        SizeOff() { }
        SizeOff(size_type s, size_type o) : size(s), offset(o) { }
        };
    size_t arrsize_ = 0;
    std::array<SizeOff,MaxNAlloc> sos_;
    T* p_ = nullptr;
    Alloc a_;
    public:

    MultAlloc() { }

    ~MultAlloc()
        {
        if(p_)
            {
            CHECK_NOT_EMPTY
            auto totsize = sos_[arrsize_-1].offset+sos_[arrsize_-1].size;
            a_.deallocate(p_,totsize);
            }
        }

    size_type
    size() const { return arrsize_; }

    size_type
    size(size_t i) const
        {
        CHECK_IND(i) 
        return sos_[i].size;
        }

    size_type constexpr
    max_nalloc() const { return MaxNAlloc; }

    void
    add(size_type size)
        {
        CHECK_NOT_FULL
        if(p_) throw std::runtime_error("Can't add to MultAlloc after allocated");
        if(arrsize_==0)
            {
            sos_[arrsize_] = SizeOff(size,0);
            }
        else
            {
            auto& prev = sos_[arrsize_-1];
            sos_[arrsize_] = SizeOff(size,prev.offset+prev.size);
            }
        ++arrsize_;
        }

    void
    allocate()
        {
        CHECK_NOT_EMPTY
        if(!p_)
            {
            auto totsize = sos_[arrsize_-1].offset+sos_[arrsize_-1].size;
            p_ = a_.allocate(totsize);
            }
        }

    pointer
    operator[](size_t i) const
        { 
        CHECK_IND(i) 
        CHECK_ALLOCATED
        CHECK_SIZE(i)
        return p_+sos_[i].offset; 
        }

    private:
    void
    check_ind(size_t i) const
        {
        if(i >= arrsize_) throw std::out_of_range("index out of range in MultAlloc");
        }
    void
    check_size(size_t i) const
        {
        if(sos_[i].size==0) throw std::out_of_range("attempted to access size zero element of MultAlloc");
        }
    void
    check_not_empty() const
        {
        if(arrsize_==0) throw std::out_of_range("MultAlloc is empty");
        }
    void
    check_not_full() const
        {
        if(arrsize_ >= max_nalloc()) throw std::out_of_range("exceeded max number in MultAlloc");
        }
    void
    check_allocated() const
        {
        if(!p_) throw std::runtime_error("MultAlloc has not been allocated");
        }
    };

#undef CHECK_IND
#undef CHECK_SIZE
#undef CHECK_NOT_FULL
#undef CHECK_NOT_EMPTY
#undef CHECK_ALLOCATED

} //namespace itensor

#endif
