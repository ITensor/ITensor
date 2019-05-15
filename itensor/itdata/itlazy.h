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
#ifndef __ITENSOR_ITLAZY_H
#define __ITENSOR_ITLAZY_H

#include "itensor/itdata/task_types.h"
#include "itensor/itdata/itdata.h"

namespace itensor {

class ITLazy
    {
    public:
    struct IndSto
        {
        IndexSet i;
        PData s;
        IndSto() { }
        IndSto(const IndexSet& i_,
               const PData& s_)
            : i(i_), s(s_) { }
        };
    using storage_type = std::vector<IndSto>;
    private:
    storage_type todo_;
    PData result_;
    public:

    ITLazy() { }

    ITLazy(const IndexSet& i1,
           const PData& s1,
           const IndexSet& i2,
           const PData& s2)
      : todo_(2) 
        { 
        todo_[0].i = i1;
        todo_[0].s = s1;
        todo_[1].i = i2;
        todo_[1].s = s2;
        }

    storage_type&
    todo() { return todo_; }

    const storage_type&
    todo() const { return todo_; }

    const IndexSet&
    iset(size_t n) const { return todo_[n].i; }

    PData&
    store(size_t n) { return todo_[n].s; }

    const PData&
    store(size_t n) const { return todo_[n].s; }

    void
    addStore(const IndexSet& is,
             const PData& pstore)
        {
        todo_.emplace_back(is,pstore);
        }

    void
    addStore(const ITLazy& other)
        {
        todo_.insert(todo_.end(),other.todo_.begin(),other.todo_.end());
        }

    bool
    hasResult() const { return static_cast<bool>(result_); }

    const PData&
    result() const { return result_; }

    void
    setResult(PData&& p) { result_ = std::move(p); }

    };

bool
hasResult(const ITLazy& Z);

PData
evaluate(ITLazy& Z);

void
doTask(Contract<Index>& C,
       ITLazy& L,
       const ITLazy& R);

void
doTask(Contract<Index>& C,
       ITLazy& L,
       const PData& R);

void
doTask(Contract<Index>& C,
       const PData& R,
       const ITLazy& L,
       ManageStore& m);

} //namespace itensor

#endif

