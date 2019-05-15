//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_MULTALLOC_H
#define __ITENSOR_MULTALLOC_H

#include <memory>
#include <vector>
#include "itensor/global.h"

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

template<typename T, size_t MaxNAlloc>
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
    size_type arrsize_ = 0;
    std::array<SizeOff,MaxNAlloc> sos_;
    std::vector<T> v_;
    //T* p_ = nullptr;
    public:

    MultAlloc() { }

    //~MultAlloc()
    //    {
    //    if(p_)
    //        {
    //        if(Global::debug1()) println("<<<< Deallocating MultAlloc, p_=",p_);
    //        CHECK_NOT_EMPTY
    //        auto totsize = sos_[arrsize_-1].offset+sos_[arrsize_-1].size;
    //        a_.deallocate(p_,totsize);
    //        }
    //    }

    size_type
    size() const { return arrsize_; }

    size_type
    data_size() const { return v_.size(); }

    size_type
    size(size_type i) const
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
        //if(p_) throw std::runtime_error("Can't add to MultAlloc after allocated");
        if(!v_.empty()) throw std::runtime_error("Can't add to MultAlloc after allocated");
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
        auto totsize = sos_[arrsize_-1].offset+sos_[arrsize_-1].size;
        v_.resize(totsize);
        }

    pointer
    operator[](size_type i)
        { 
        CHECK_IND(i) 
        CHECK_ALLOCATED
        CHECK_SIZE(i)
#ifdef DEBUG
        if(sos_[i].offset+sos_[i].size > v_.size()) throw std::out_of_range("data out of range in MultAlloc");
#endif
        return v_.data()+sos_[i].offset; 
        }

    private:
    void
    check_ind(size_type i) const
        {
        if(i >= arrsize_) throw std::out_of_range("index out of range in MultAlloc");
        }
    void
    check_size(size_type i) const
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
        if(v_.empty()) throw std::runtime_error("MultAlloc has not been allocated");
        }
    };

#undef CHECK_IND
#undef CHECK_SIZE
#undef CHECK_NOT_FULL
#undef CHECK_NOT_EMPTY
#undef CHECK_ALLOCATED

} //namespace itensor

#endif
