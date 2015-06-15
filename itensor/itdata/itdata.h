//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_H
#define __ITENSOR_ITDATA_H

#include "itensor/types.h"
#include "itensor/util/error.h"
#include "itensor/itdata/storage_types.h"

namespace itensor {

struct ITData;
using PData = std::shared_ptr<ITData>;
using CPData = std::shared_ptr<const ITData>;

template<typename List>
struct FuncBaseT : FuncBaseT<typename List::Next>
    {
    using T = typename List::Type;
    using FuncBaseT<typename List::Next>::applyTo;

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
    public:

    T d;

    template<typename... VArgs>
    ITWrap(VArgs&&... vargs) : d(std::forward<VArgs>(vargs)...) 
        { 
        static_assert(containsType<StorageTypes,T>{},"Data type not in list of registered storage types");
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
    PData *parg1_ = nullptr;
    const CPData *parg2_ = nullptr;
    const ITData *arg2_ = nullptr;
    Action action_ = None;
    PData nd_;

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
                //*pdata_ = std::make_shared<T>(*olda1);
                *pdata_ = std::make_shared<ITWrap<T>>(*olda1);
                }
            return *(static_cast<T*>(pdata_->get()));
            }
        };

    public:

    ManageStore() { }

    ManageStore(PData *parg1)
        : parg1_(parg1)
        { }

    ManageStore(PData *parg1, const ITData *arg2)
        : parg1_(parg1), arg2_(arg2)
        { }

    ManageStore(PData *parg1, const CPData *parg2)
        : parg1_(parg1), parg2_(parg2), arg2_(parg2->get())
        { }

    ManageStore(ManageStore&& o)
        :
        parg1_(o.parg1_),
        parg2_(o.parg2_),
        arg2_(o.arg2_),
        action_(o.action_),
        nd_(std::move(o.nd_))
        { 
        o.parg1_ = nullptr;
        o.parg2_ = nullptr;
        o.arg2_ = nullptr;
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
    hasPArg1() const { return bool(parg1_); }

    bool
    hasPArg2() const { return bool(parg2_); }

    bool
    hasArg2() const { return bool(arg2_); }

    PData&
    parg1() 
        { 
#ifdef DEBUG
        if(!parg1_) Error("Attempt to dereference nullptr");
#endif
        return *parg1_; 
        }


    const CPData&
    parg2() 
        { 
#ifdef DEBUG
        if(!parg2_) Error("Attempt to dereference nullptr");
#endif
        return *parg2_; 
        }

    const ITData&
    arg2() 
        { 
#ifdef DEBUG
        if(!arg2_) Error("Attempt to dereference nullptr");
#endif
        return *arg2_; 
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

    private:

    void
    updateArg1();

    };

template <typename StorageT, typename... VArgs>
StorageT* ManageStore::
makeNewData(VArgs&&... vargs)
    {
    if(!parg1_) Error("Can't call makeNewData with const-only access to first arg");
    action_ = AssignNewData;
    auto newdat = std::make_shared<ITWrap<StorageT>>(std::forward<VArgs>(vargs)...);
    auto* ret = newdat.get();
    nd_ = std::move(newdat);
    return &(ret->d);
    }

void inline ManageStore::
assignPointerRtoL() 
    { 
    if(!parg2_) Error("No second pointer provided for action AssignPointerRtoL");
    action_ = AssignPointerRtoL; 
    }

void inline ManageStore::
updateArg1()
    {
    if(!parg1_) return;
    //println("In updateArg1, arg1_ points to ",arg1_->get());
    if(action_ == AssignNewData)
        {
        //println("Doing AssignNewData");
        *parg1_ = std::move(nd_);
        }
    else if(action_ == AssignPointerRtoL)
        {
        //println("Doing AssignPointerRtoL");
        *parg1_ = std::const_pointer_cast<ITData,const ITData>(*parg2_);
        }
    }

ManageStore::UniqueRef inline ManageStore::
modifyData()
    {
    if(!parg1_) Error("Can't modify const data");
    return UniqueRef(parg1_);
    }

template<typename T>
T* ManageStore::
modifyData(const T& d)
    {
    if(!parg1_) Error("Can't modify const data");
    if(!(parg1_->unique())) 
        {
        auto* olda1 = static_cast<ITWrap<T>*>(parg1_->get());
        *parg1_ = std::make_shared<ITWrap<T>>(olda1->d);
        }
    auto* a1 = static_cast<ITWrap<T>*>(parg1_->get());
    return &(a1->d);
    }


} // namespace itensor

#endif
