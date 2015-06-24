//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITLAZY_H
#define __ITENSOR_ITLAZY_H

#include "itensor/itdata/task_types.h"
#include "itensor/itdata/itdata.h"

namespace itensor {

class ITReal;
class ITCplx;

class ITLazy
    {
    public:
    struct IndSto
        {
        IndexSet i;
        PData s;
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

    PData&
    result() { return result_; }

    void
    addStore(const IndexSet& is,
             const PData& pstore)
        {
        todo_.emplace_back(is,pstore);
        }

    };

void inline
evaluate(ITLazy& L, ManageStore& m)
    {
    //0. Check if already evaluated
    if(L.result()) 
        {
        m.pointTo(L.result());
        return;
        }
    //1. Reorder intermediate indices
    //2. Do chain of contractions
    }

void inline
doTask(Contract<Index>& C,
       ITLazy& L,
       const ITReal& R,
       ManageStore& m)
    {
    contractIS(C.Lis,C.Ris,C.Nis);
    L.addStore(C.Nis,m.parg2());
    }

void inline
doTask(Contract<Index>& C,
       const ITReal& R,
       const ITLazy& L,
       ManageStore& m)
    {
    //TODO: create new index set
    auto* pn = makeNewData(L);
    Error("Need access to arg1 pointer");
    //TODO: need to be able to access pointer
    //      to arg1 here
    //pn->addStore(new_indexset,m.parg1());
    }

} //namespace itensor

#endif

