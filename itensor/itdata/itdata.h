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
#ifndef __ITENSOR_ITDATA_H
#define __ITENSOR_ITDATA_H

#include <memory>
#include "itensor/types.h"
#include "itensor/util/error.h"
#include "itensor/util/timers.h"
#include "itensor/itdata/storage_types.h"

namespace itensor {

struct ITData;

using PData = std::shared_ptr<ITData>;

struct CPData  //logically const ITData smart pointer
    {
    PData& p;

    CPData(PData& p_) : p(p_) { }

    explicit operator bool() const { return static_cast<bool>(p); }

    ITData&
    operator*();

    const ITData&
    operator*() const;

    ITData*
    operator->();

    const ITData*
    operator->() const;
    };

template<typename TList>
struct FuncBaseT : FuncBaseT<popFront<TList>>
    {
    using T = frontType<TList>;
    using FuncBaseT<popFront<TList>>::applyTo;

    void virtual
    applyTo(const T& t) = 0;

    void virtual
    applyTo(T& t) = 0;
    };
template<>
struct FuncBaseT<TypeList<>>
    {
    void
    applyTo() { }
    };

using FuncBase = FuncBaseT<StorageTypes>;


struct ITData
    {
    ITData() { }
    virtual ~ITData() { }

    PData virtual
    clone() const = 0;

    void virtual
    plugInto(FuncBase& f) const = 0;

    void virtual
    plugInto(FuncBase& f) = 0;
    };

template<typename T>
class ITWrap : public ITData
    {
    static_assert(containsType<StorageTypes,stdx::decay_t<T>>{},"Data type not in list of registered storage types");
    public:

    T d;

    template<typename... VArgs>
    ITWrap(VArgs&&... vargs) : d(std::forward<VArgs>(vargs)...) 
        { 
        }

    virtual ~ITWrap() { }

    private:
    
    PData
    clone() const final 
        { 
        return std::make_shared<ITWrap<T>>(d);
        }

    void
    plugInto(FuncBase& f) const final
        {
        f.applyTo(d);
        }
        
    void
    plugInto(FuncBase& f) final
        {
        f.applyTo(d);
        }
    };


class ManageStore
    {
    enum Action
        {
        None,
        AssignNewData,
        AssignPointerRtoL
        };

    PData *pparg1_ = nullptr;
    PData *pparg2_ = nullptr;
    Action action_ = None;
    PData nd_;

    class UniqueRef;

    public:

    ManageStore() { }

    ManageStore(PData *pparg1)
      : pparg1_(pparg1)
        { }

    ManageStore(PData *pparg1, PData *pparg2)
      : pparg1_(pparg1), pparg2_(pparg2)
        { }

    ManageStore(ManageStore&& o)
      : pparg1_(o.pparg1_),
        pparg2_(o.pparg2_),
        action_(o.action_),
        nd_(std::move(o.nd_))
        { 
        o.pparg1_ = nullptr;
        o.pparg2_ = nullptr;
        o.action_ = None;
        }

    ~ManageStore()
        {
        updateArg1();
        }

    ManageStore(const ManageStore&) = delete;

    ManageStore&
    operator=(const ManageStore&) = delete;

    bool
    hasPArg1() const { return bool(pparg1_); }

    bool
    hasPArg2() const { return bool(pparg2_); }

    void
    setparg1(PData *pparg1)
        {
        pparg1_ = pparg1;
        }

    void
    setparg2(PData *pparg2)
        {
        pparg2_ = pparg2;
        }

    PData&
    parg1() 
        { 
#ifdef DEBUG
        if(!pparg1_) Error("Attempt to dereference nullptr");
#endif
        return *pparg1_; 
        }

    PData&
    parg2() 
        { 
#ifdef DEBUG
        if(!pparg2_) Error("Attempt to dereference nullptr");
#endif
        return *pparg2_; 
        }

    //Can be used as StorageType& sref = mp.modifyData();
    //The RefHelper return type will deduce StorageType
    //and convert to a StorageType&
    UniqueRef
    modifyData();

    //This returns a pointer because if it returned a reference
    //it is too easy to copy the data type by writing
    //auto nd = modifyData(...);
    template<typename T>
    T*
    modifyData(const T& d);

    //This returns a pointer because if it returned a reference
    //it is too easy to copy the data type by writing
    //auto nd = makeNewData(...);
    template <typename StorageT, typename... Args>
    StorageT*
    makeNewData(Args&&... args);

    PData&
    newData() { return nd_; }

    void
    assignPointerRtoL();

    void
    pointTo(const PData& p);

    private:

    void
    updateArg1();

    class UniqueRef
        {
        PData* pdata_ = nullptr;
        public:

        UniqueRef(PData* pdata) : pdata_(pdata) { }

        template<typename T>
        operator T&()
            {
            if(!(pdata_->unique())) 
                {
                auto* olda1 = static_cast<T*>(pdata_->get());
                *pdata_ = std::make_shared<ITWrap<T>>(*olda1);
                }
            return *(static_cast<T*>(pdata_->get()));
            }
        };
    };

template <typename StorageT, typename... VArgs>
StorageT* ManageStore::
makeNewData(VArgs&&... vargs)
    {
    //if(!pparg1_) Error("Can't call makeNewData with const-only access to first arg");
    action_ = AssignNewData;
    auto newdat = std::make_shared<ITWrap<StorageT>>(std::forward<VArgs>(vargs)...);
    auto* ret = newdat.get();
    nd_ = std::move(newdat);
    return &(ret->d);
    }

void inline ManageStore::
assignPointerRtoL() 
    { 
    //if(!pparg2_) Error("No second pointer provided for action AssignPointerRtoL");
    action_ = AssignPointerRtoL; 
    }

void inline ManageStore::
pointTo(const PData& p) 
    { 
    action_ = AssignNewData; 
    nd_ = p;
    }

void inline ManageStore::
updateArg1()
    {
    if(!pparg1_) return;
    //println("In updateArg1, arg1_ points to ",arg1_->get());
    if(action_ == AssignNewData)
        {
        //println("Doing AssignNewData");
        *pparg1_ = std::move(nd_);
        }
    else if(action_ == AssignPointerRtoL)
        {
        //println("Doing AssignPointerRtoL");
        //*pparg1_ = std::const_pointer_cast<ITData,const ITData>(*pparg2_);
        *pparg1_ = *pparg2_;
        }
    }

ManageStore::UniqueRef inline ManageStore::
modifyData()
    {
    //if(!pparg1_) Error("Can't modify const data");
    return UniqueRef(pparg1_);
    }

template<typename T>
T* ManageStore::
modifyData(const T& d)
    {
    //if(!pparg1_) Error("Can't modify const data");
    if(!(pparg1_->unique())) 
        {
        auto* olda1 = static_cast<ITWrap<T>*>(pparg1_->get());
        *pparg1_ = std::make_shared<ITWrap<T>>(olda1->d);
        }
    auto* a1 = static_cast<ITWrap<T>*>(pparg1_->get());
    return &(a1->d);
    }


inline ITData& CPData::
operator*() { return *p; }

inline const ITData& CPData::
operator*() const { return *p; }

inline ITData* CPData::
operator->() { return p.get(); }

inline const ITData* CPData::
operator->() const { return p.get(); }

template<typename T>
const char*
typeNameOf(T const& t) { return "[unknown]"; }


} // namespace itensor

#endif
