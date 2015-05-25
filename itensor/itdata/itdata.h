//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_H
#define __ITENSOR_ITDATA_H
#include "itensor/global.h"
#include "itensor/detail/algs.h"
#include "itensor/itdata/storage_types.h"
//#include "itensor/detail/call_rewrite.h"

#define LPAREN (
#define RPAREN )

namespace itensor {

struct ITData;
using PData = std::shared_ptr<ITData>;
using CPData = std::shared_ptr<const ITData>;

struct FuncBase
    {
    FuncBase() { }
    virtual ~FuncBase() { }

    REGISTER_TYPES(void virtual applyTo LPAREN, &d RPAREN =0;)
    REGISTER_TYPES(void virtual applyTo LPAREN, const& d RPAREN =0;)

    template <typename T>
    void
    applyTo(T& t) { throw ITError("ITData subtype not registered."); }
    };

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

template <class DType>
struct ITDataType : ITData
    {
    DType d;

    ITDataType() { }

    template<typename... VArgs>
    ITDataType(VArgs&&... vargs)
        : d(std::forward<VArgs>(vargs)...)
        { }

    virtual ~ITDataType() { }

    private:
    
    PData
    clone() const final 
        { 
        return std::make_shared<ITDataType<DType>>(d);
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

//////////////////
//////////////////

struct Void { };

class ManagePtr
    {
    enum Action
        {
        None,
        AssignNewData,
        AssignPointerRtoL
        };
    PData *parg1_ = nullptr;
    const PData *parg2_ = nullptr;
    const ITData *arg2_ = nullptr;
    Action action_ = None;
    PData nd_;
    public:

    ManagePtr() { }

    ManagePtr(PData *parg1)
        : parg1_(parg1)
        { }

    ManagePtr(PData *parg1, const ITData *arg2)
        : parg1_(parg1), arg2_(arg2)
        { }

    ManagePtr(PData *parg1, const PData *parg2)
        : parg1_(parg1), parg2_(parg2), arg2_(parg2->get())
        { }

    ManagePtr(ManagePtr&& o)
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

    ~ManagePtr()
        {
        updateArg1();
        }

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
        if(!parg1_) throw std::runtime_error("Attempt to dereference nullptr");
#endif
        return *parg1_; 
        }


    const PData&
    parg2() 
        { 
#ifdef DEBUG
        if(!parg2_) throw std::runtime_error("Attempt to dereference nullptr");
#endif
        return *parg2_; 
        }

    const ITData&
    arg2() 
        { 
#ifdef DEBUG
        if(!arg2_) throw std::runtime_error("Attempt to dereference nullptr");
#endif
        return *arg2_; 
        }


    //This returns a pointer because otherwise
    //it is too easy to copy the data type by writing
    //auto nd = modifyData(...); (should be auto& if reference returned)
    template <typename T>
    T*
    modifyData(const T& d);

    //This returns a pointer because otherwise
    //it is too easy to copy the data type by writing
    //auto nd = makeNewData(...); (should be auto& if reference returned)
    template <typename StorageT, typename... Args>
    StorageT*
    makeNewData(Args&&... args);

    //template <typename StorageT>
    //void
    //setNewData(std::shared_ptr<ITDataType>&& nd);

    PData&
    newData() { return nd_; }

    void
    assignPointerRtoL();

    private:

    void
    updateArg1();

    };

//OneArg and TwoArgs are "policy classes" 
//for customizing implementation of RegisterTask
struct OneArg
    {
    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, const D& d, ManagePtr& mp);

    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, D& d, ManagePtr& mp);
    };

struct TwoArgs
    {
    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, const D& d, ManagePtr& mp);

    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, D& d, ManagePtr& mp);
    };

template <typename NArgs, typename Task, typename Return>
class RegisterTask : public FuncBase
    {
    Task t_;
    ManagePtr mp_;
    Return ret_;
    public:
    using return_type = Return;

    RegisterTask(Task&& t) :
        t_(std::move(t))
        { }

    template<typename... MPArgs>
    RegisterTask(Task&& t,
                 MPArgs&&... mpargs) :
        t_(std::move(t)),
        mp_(std::forward<MPArgs>(mpargs)...)
        { }

    RegisterTask(RegisterTask&& o) :
        t_(std::move(o.t_)), 
        mp_(std::move(o.mp_)),
        ret_(std::move(o.ret_))
        { }

    virtual ~RegisterTask() { }

    operator return_type() { return std::move(ret_); }

    REGISTER_TYPES(void applyTo LPAREN, &d RPAREN final { ret_ = std::move(NArgs().template call<return_type>(t_,d,mp_)); } )
    REGISTER_TYPES(void applyTo LPAREN, const&d RPAREN final { ret_ = std::move(NArgs().template call<return_type>(t_,d,mp_)); } )
    };

template <typename Task, typename D1, typename Return>
struct CallWrap : FuncBase
    {
    CallWrap(Task& t, D1& arg1, ManagePtr& mp) 
        : t_(t), arg1_(arg1), parg1_(nullptr), mp_(mp) { }
    CallWrap(Task& t, D1& arg1, PData& parg1, ManagePtr& mp) 
        : t_(t), arg1_(arg1), parg1_(&parg1), mp_(mp) { }

    Return
    getReturn() { return std::move(ret_); }

    private:

    template<typename D2>
    void
    applyToImpl(const D2& d2);
    
    Task& t_;
    D1& arg1_;
    PData* parg1_;
    Return ret_;
    ManagePtr& mp_;

    public:
    REGISTER_TYPES(void applyTo LPAREN, &d RPAREN final { applyToImpl(d); })
    REGISTER_TYPES(void applyTo LPAREN, const&d RPAREN final { applyToImpl(d); })
    };



//
// Implementations
//

namespace detail {

//If ActualRet!=void this gets called
template<typename Ret, typename ActualRet, typename... VArgs>
struct FixRet
    {
    Ret
    operator()(VArgs&&... vargs)
        {
        static_assert(!std::is_same<Ret,Void>::value,
                      "No return type specified in call to doTask");
        static_assert(std::is_same<Ret,ActualRet>::value,
                      "Mismatched return type specified in call to doTask");
        return doTask(std::forward<VArgs>(vargs)...);
        }
    };
//If ActualRet==void this gets called
template<typename Ret, typename... VArgs>
struct FixRet<Ret,void,VArgs...>
    {
    Ret
    operator()(VArgs&&... vargs)
        {
        doTask(std::forward<VArgs>(vargs)...);
        return Ret{};
        }
    };

//(4) Least preferred version which always compiles
template <typename Ret, typename Task, typename Storage>
Ret
callDoTask_Impl(long l, Task& t, Storage& s)
    {
    throw std::runtime_error("1 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
//(3) This is preferred more than version (4) above
// since 0->int requires no cast.
// Template substitution will fail (not a compile-time error)
// if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret, typename Task, typename Storage>
auto 
callDoTask_Impl(int i, Task& t, Storage& s)
    -> std::conditional_t<std::is_same<decltype(doTask(t,s)),void>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,s));
    return FixRet<Ret,ActualRet,Task&,Storage&>()(t,s);
    }
//(2) Less preferred version which always compiles and passes
//the call on to callDoTask_Impl without the ManagePtr& argument
template <typename Ret, typename Task, typename Storage>
Ret
callDoTask_Impl(long l, Task& t, Storage& s, ManagePtr& mp)
    {
    return callDoTask_Impl<Ret>(0,t,s);
    }
//(1) This is the preferred version since 0->int requires no cast
//Template substitution will fail (not a compile-time error)
//if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret, typename Task, typename Storage>
auto 
callDoTask_Impl(int i, Task& t, Storage& s, ManagePtr& mp)
    -> std::conditional_t<std::is_same<decltype(doTask(t,s,mp)),void>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,s,mp));
    return FixRet<Ret,ActualRet,Task&,Storage&,ManagePtr&>()(t,s,mp);
    }
//This version of callDoTask attempts "return doTask(t,s);"
//- If doTask(t,s) not defined, throws an exception
//  if called (this converts what would otherwise be 
//  a compile-time error into a run-time error)
//- If doTask(t,s) defined but returns void, 
//  returns Ret{} instead
template<typename Ret, typename Task, typename Storage>
Ret
callDoTask(Task& t, Storage& s)
    {
    //Skip straight to (3) above
    return callDoTask_Impl<Ret,Task,Storage>(0,t,s);
    }
//This version of callDoTask attempts "return doTask(t,s,mp);"
//- If doTask(t,s,mp) not defined, tries calling doTask(t,s)
//  then throws an exception if that is not defined
//  (this converts what would otherwise be 
//  a compile-time error into a run-time error)
//- If doTask(t,s,mp) or doTask(t,s) defined but 
//  either returns void, returns Ret{} instead
template<typename Ret, typename Task, typename Storage>
Ret
callDoTask(Task& t, Storage& s, ManagePtr& mp)
    {
    //First try calling function labeled (1) above
    return callDoTask_Impl<Ret,Task,Storage>(0,t,s,mp);
    }

/////////////////////////////////////////////////////

//(4) Least preferred version calls doTask(Task,D&)
//and makes sure pdat is unique first ("copy on write")
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Impl(Task& t, const D& cd, PData& pdat,long)
    {
    //println("--> Calling solo (1 param, no ManagePtr)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pd = static_cast<ITDataType<D>*>(pdat.get());
    return callDoTask<Ret>(t,pd->d);
    }
//(3) Preferred version since 0->int requires no conversion
//Attempts to call doTask(Task,const D&) if defined
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Impl(Task& t, const D& cd, PData& pdat,int)
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T&) to get substitution failure (SFINAE) if no such call exists.
    -> std::conditional_t<std::is_same<decltype(doTask(t,cd)),void>::value,Ret,Ret>
    {
    //println("--> Not calling solo (1 param, no ManagePtr)");
    return callDoTask<Ret>(t,cd);
    }
//(2) Less preferred version calls one of the cloneDoTask_Impl
//versions (which one depends on further tests - see above)
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_MPImpl(Task& t, const D& cd, D& d, ManagePtr& mp,long)
    {
    //Give up on trying to pass ManagePtr argument here - 
    //this assumes no need to have ManagePtr if definitely
    //going to modify storage in-place
    return cloneDoTask_Impl<Ret>(t,d,mp.parg1(),0);
    }
//(1) Preferred version since 0->int requires no conversion
//Attempts to call doTask(Task,const D&, ManagePtr&) if defined
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_MPImpl(Task& t, const D& cd, D& d, ManagePtr& mp,int)
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T&) to get substitution failure (SFINAE) if no such call exists.
    -> std::conditional_t<std::is_same<decltype(doTask(t,cd,mp)),void>::value,Ret,Ret>
    {
    //println("--> Not calling solo (1 param + optional ManagePtr)");
    return callDoTask<Ret>(t,cd,mp);
    }
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask(Task& t, D& d, ManagePtr& mp)
    {
    return cloneDoTask_MPImpl<Ret>(t,d,d,mp,0);
    }

/////////////////////

//(4) Least preferred version which always compiles
template <typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask_Impl(long l, Task& t, D1& d1, const D2& d2)
    {
    throw std::runtime_error("2 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
//(3) This is preferred more than version (4) above
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(int i, Task& t, D1& d1, const D2& d2)
    -> std::conditional_t<std::is_same<decltype(doTask(t,d1,d2)),void>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,d1,d2));
    return FixRet<Ret,ActualRet,Task&,D1&,const D2&>()(t,d1,d2);
    }
//(2) Less preferred version which always compiles and passes
//the call on to callDoTask_Impl (3) without the ManagePtr& argument
template <typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask_Impl(long l, Task& t, D1& d1, const D2& d2, ManagePtr& mp)
    {
    return callDoTask_Impl<Ret>(0,t,d1,d2);
    }
//(1)
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(int i, Task& t, D1& d1, const D2& d2, ManagePtr& mp)
    -> std::conditional_t<std::is_same<decltype(doTask(t,d1,d2,mp)),void>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,d1,d2,mp));
    return FixRet<Ret,ActualRet,Task&,D1&,const D2&,ManagePtr&>()(t,d1,d2,mp);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask(Task& t, D1& d1, const D2& d2)
    {
    //Skip straight to (3) above
    return callDoTask_Impl<Ret,Task,D1,D2>(0,t,d1,d2);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask(Task& t, D1& d1, const D2& d2, ManagePtr& mp)
    {
    //First try calling function labeled (1) above
    return callDoTask_Impl<Ret,Task,D1,D2>(0,t,d1,d2,mp);
    }

/////////////////////

//(4)
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_Impl(Task& t, const D1& cd1, D1& d1, const D2& d2, PData& p1,long)
    {
    //println("--> Calling solo (2 params)");
    if(!p1.unique()) p1 = p1->clone();
    auto* pd1 = static_cast<ITDataType<D1>*>(p1.get());
    return callDoTask<Ret>(t,pd1->d,d2);
    }
//(3)
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_Impl(Task& t, const D1& cd1, D1& d1, const D2& d2, PData& p1,int) 
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T1&,const T2&) to get substitution failure (SFINAE) 
    //if no such call exists.
    -> std::conditional_t<std::is_same<decltype(doTask(t,cd1,d2)),void>::value,Ret,Ret>
    {
    //println("--> Not calling solo (2 param, no ManagePtr)");
    return callDoTask<Ret>(t,cd1,d2);
    }
//(2)
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_MPImpl(Task& t, const D1& cd1, D1& d1, const D2& d2, ManagePtr& mp,long)
    {
    return cloneDoTask_Impl<Ret,Task,D1,D2>(t,cd1,d1,d2,mp.parg1(),0);
    }
//(1)
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_MPImpl(Task& t, const D1& cd1, D1& d1, const D2& d2, ManagePtr& mp,int) 
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T1&,const T2&) to get substitution failure (SFINAE) 
    //if no such call exists.
    -> std::conditional_t<std::is_same<decltype(doTask(t,cd1,d2,mp)),void>::value,Ret,Ret>
    {
    //println("--> Not calling solo (2 param + optional ManagePtr)");
    return callDoTask<Ret>(t,cd1,d2,mp);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask(Task& t, D1& d1, const D2& d2, ManagePtr& mp)
    {
    return cloneDoTask_MPImpl<Ret,Task,D1,D2>(t,d1,d1,d2,mp,0);
    }

} //namespace detail

template <typename StorageT, typename... VArgs>
StorageT* ManagePtr::
makeNewData(VArgs&&... vargs)
    {
    if(!parg1_) Error("Can't call makeNewData with const-only access to first arg");
    action_ = AssignNewData;
    auto newdat = std::make_shared<ITDataType<StorageT>>(std::forward<VArgs>(vargs)...);
    auto* ret = newdat.get();
    nd_ = std::move(newdat);
    return &(ret->d);
    }

//template <typename ITDataType>
//void ManagePtr::
//setNewData(std::shared_ptr<ITDataType>&& nd)
//    {
//    nd_ = std::move(nd);
//    action_ = AssignNewData;
//    }

void inline ManagePtr::
assignPointerRtoL() 
    { 
    if(!parg2_) Error("No second pointer provided for action AssignPointerRtoL");
    action_ = AssignPointerRtoL; 
    }

void inline ManagePtr::
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
        *parg1_ = *parg2_;
        }
    }

template<typename T>
T* ManagePtr::
modifyData(const T& d)
    {
    if(!parg1_) Error("Can't modify const data");
    if(!(parg1_->unique())) 
        {
        *parg1_ = (*parg1_)->clone();
        }
    auto* pa1 = static_cast<ITDataType<T>*>(parg1_->get());
    return &(pa1->d);
    }

template<typename Ret, typename Task, typename D>
Ret OneArg::
call(Task& t, const D& d, ManagePtr& mp)
    {
    return detail::callDoTask<Ret>(t,d,mp);
    }

template<typename Ret, typename Task, typename D>
Ret OneArg::
call(Task& t, D& d, ManagePtr& mp)
    {
    return detail::cloneDoTask<Ret>(t,d,mp);
    }

template<typename Ret, typename Task, typename D>
Ret TwoArgs::
call(Task& t, const D& d, ManagePtr& mp)
    {
    assert(mp.hasArg2());
    CallWrap<Task,const D,Ret> w(t,d,mp);
    mp.arg2().plugInto(w);
    return w.getReturn();
    }

template<typename Ret, typename Task, typename D>
Ret TwoArgs::
call(Task& t, D& d, ManagePtr& mp)
    {
    assert(mp.hasPArg1());
    assert(mp.hasArg2());
    CallWrap<Task,D,Ret> w(t,d,mp.parg1(),mp);
    mp.arg2().plugInto(w);
    return w.getReturn();
    }

template <typename Task, typename D1, typename Return>
template<typename D2>
void CallWrap<Task,D1,Return>::
applyToImpl(const D2& d2)
    { 
    if(parg1_) ret_ = detail::cloneDoTask<Return>(t_,arg1_,d2,mp_); 
    else       ret_ = detail::callDoTask<Return>(t_,arg1_,d2,mp_); 
    }


//////
////// doTask methods
//////


template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       const ITData& arg)
    {
    RegisterTask<OneArg,Task,ReturnType> r(std::move(t));
    arg.plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       const ITData& arg)
    {
    doTask<Void,Task>(std::forward<Task>(t),arg);
    return t;
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       const CPData& arg)
    {
    RegisterTask<OneArg,Task,ReturnType> r(std::move(t));
    arg->plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       const CPData& arg)
    {
    doTask<Void,Task>(std::forward<Task>(t),arg);
    return t;
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       PData& arg)
    {
    RegisterTask<OneArg,Task,ReturnType> r(std::move(t),&arg);
    arg->plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       PData& arg)
    {
    doTask<Void,Task>(std::forward<Task>(t),arg);
    return t;
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       const PData& arg1,
       const PData& arg2)
    {
    RegisterTask<TwoArgs,Task,ReturnType> r(std::move(t),&arg1,&arg2);
    arg1->plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       const PData& arg1,
       const PData& arg2)
    {
    doTask<Void,Task>(std::forward<Task>(t),arg1,arg2);
    return t;
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       PData& arg1,
       const PData& arg2)
    {
    RegisterTask<TwoArgs,Task,ReturnType> r(std::move(t),&arg1,&arg2);
    arg1->plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       PData& arg1,
       const PData& arg2)
    {
    doTask<Void,Task>(std::forward<Task>(t),arg1,arg2);
    return t;
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       PData& arg1,
       const ITData& arg2)
    {
    RegisterTask<TwoArgs,Task,ReturnType> r(std::move(t),&arg1,&arg2);
    arg1->plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       PData& arg1,
       const ITData& arg2)
    {
    doTask<Void,Task>(std::forward<Task>(t),arg1,arg2);
    return t;
    }

} //namespace itensor

#undef LPAREN
#undef RPAREN
#undef REGISTER_TYPES

#endif

