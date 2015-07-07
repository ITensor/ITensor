//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITLAZY_H
#define __ITENSOR_ITLAZY_H

#include "itensor/itdata/task_types.h"
#include "itensor/itdata/itdata.h"
#include "itensor/itdata/synchronized.h"

namespace itensor {

class ITLazy
    {
    public:
    struct IndSto
        {
        IndexSet i;
        CPData s;
        IndSto() { }
        IndSto(const IndexSet& i_,
               const CPData& s_)
            : i(i_), s(s_) { }
        };
    using storage_type = std::vector<IndSto>;
    private:
    storage_type todo_;
    mutable Synchronized<PData> result_;
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

    void
    addStore(const IndexSet& is,
             const CPData& pstore)
        {
        todo_.emplace_back(is,pstore);
        }

    bool
    hasResult() const { return static_cast<bool>(result_); }

    PData
    result() const { return result_.get(); }

    void
    setResult(PData&& p) const { result_.set(std::move(p)); }

    };

PData inline
evaluate(const ITLazy& L)
    {
    //0. Check if already evaluated
    if(L.hasResult()) return L.result();
    //1. Reorder intermediate indices
    //2. Do chain of contractions
    auto p = std::make_shared<ITWrap<ITReal>>(10,2);
    L.setResult(std::move(p));
    return L.result();
    }

void inline
doTask(Contract<Index>& C,
       ITLazy& L,
       const ITLazy& R)
    {
    }

void inline
doTask(Contract<Index>& C,
       ITLazy& L,
       const CPData& R,
       ManageStore& m)
    {
    contractIS(C.Lis,C.Ris,C.Nis);
    L.addStore(C.Nis,m.parg2());
    }


void inline
doTask(Contract<Index>& C,
       const CPData& R,
       const ITLazy& L,
       ManageStore& m)
    {
    //TODO: create new index set
    //auto* pn = m.makeNewData(L);
    Error("Need access to arg1 pointer");
    //TODO: need to be able to access pointer
    //      to arg1 here
    //pn->addStore(new_indexset,m.parg1());
    }

} //namespace itensor

#endif

