//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_H
#define __ITENSOR_ITDATA_H
#include "itensor/global.h"
#include "itensor/detail/algs.h"
#include "itensor/itdata/storage_types.h"

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
    applyTo(T& t) { Error("ITData subtype not registered."); }
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

class ManagePtr
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
    public:

    ManagePtr() { }

    ManagePtr(PData *parg1)
        : parg1_(parg1)
        { }

    ManagePtr(PData *parg1, const ITData *arg2)
        : parg1_(parg1), arg2_(arg2)
        { }

    ManagePtr(PData *parg1, const CPData *parg2)
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
        *parg1_ = std::const_pointer_cast<ITData,const ITData>(*parg2_);
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

//////////////////

struct Void { };


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
    public:
    using task_type = std::remove_reference_t<Task>;
    using return_type = std::remove_reference_t<Return>;
    private:
    task_type task_;
    ManagePtr mp_;
    return_type ret_;
    public:

    RegisterTask(task_type&& t) :
        task_(std::move(t))
        { }

    template<typename... MPArgs>
    RegisterTask(task_type&& t,
                 MPArgs&&... mpargs) :
        task_(std::move(t)),
        mp_(std::forward<MPArgs>(mpargs)...)
        { }

    RegisterTask(RegisterTask&& o) :
        task_(std::move(o.task_)), 
        mp_(std::move(o.mp_)),
        ret_(std::move(o.ret_))
        { }

    virtual ~RegisterTask() { }

    return_type&
    getReturn() { return ret_; }

    task_type&
    getTask() { return task_; }

    private:

    REGISTER_TYPES(void applyTo LPAREN, &d RPAREN final { applyToImpl(d); } )
    REGISTER_TYPES(void applyTo LPAREN, const&d RPAREN final { applyToImpl(d); } )

    template<typename D>
    void
    applyToImpl(const D& d)
        {
        ret_ = std::move(NArgs().template call<return_type>(task_,d,mp_));
        }
    template<typename D>
    void
    applyToImpl(D& d)
        {
        ret_ = std::move(NArgs().template call<return_type>(task_,d,mp_));
        }
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
    Error("1 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
//(3) This is preferred more than version (4) above
// since 0->int requires no cast.
// Template substitution will fail (not a compile-time error)
// if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret, typename Task, typename Storage>
auto 
callDoTask_Impl(int i, Task& t, Storage& s)
    -> std::conditional_t<std::is_void<decltype(doTask(t,s))>::value,Ret,Ret>
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
    -> std::conditional_t<std::is_void<decltype(doTask(t,s,mp))>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,s,mp));
    return FixRet<Ret,decltype(doTask(t,s,mp)),Task&,Storage&,ManagePtr&>()(t,s,mp);
    }
//This version of callDoTask attempts "return doTask(t,s,mp);"
//- If doTask(t,s,mp) not defined, tries calling doTask(t,s)
//  then aborts if that is not defined
//  (this converts what would otherwise be 
//  a compile-time error into a run-time error)
//- If doTask(t,s,mp) or doTask(t,s) defined but 
//  either returns void, returns Ret{} instead
template<typename Ret, typename Task, typename Storage>
Ret
callDoTask(Task& t, Storage& s, ManagePtr& mp)
    {
    return callDoTask_Impl<Ret,Task,Storage>(0,t,s,mp);
    }

/////////////////////////////////////////////////////

//(2-fail)
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Case2ConstNoMP(Task& t, const D& cd, ManagePtr& mp,long)
    {
    if(!mp.parg1().unique()) mp.parg1() = mp.parg1()->clone();
    auto* pd = static_cast<ITDataType<D>*>(mp.parg1().get());
    return callDoTask<Ret>(t,pd->d,mp);
    }
//(2-success)
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Case2ConstNoMP(Task& t, const D& cd, ManagePtr& mp,int)
    -> std::conditional_t<std::is_void<decltype(doTask(t,cd))>::value,Ret,Ret>
    {
    return FixRet<Ret,decltype(doTask(t,cd)),Task&,const D&>()(t,cd);
    }
//(1-fail)
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Case1ConstMP(Task& t, const D& cd, ManagePtr& mp,long)
    {
    return cloneDoTask_Case2ConstNoMP<Ret>(t,cd,mp,0);
    }
//(1-success) Preferred version since 0->int requires no conversion
//Attempts to call doTask(Task,const D&, ManagePtr&) if defined
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Case1ConstMP(Task& t, const D& cd, ManagePtr& mp,int)
    -> std::conditional_t<std::is_void<decltype(doTask(t,cd,mp))>::value,Ret,Ret>
    {
    return FixRet<Ret,decltype(doTask(t,cd,mp)),Task&,const D&,ManagePtr&>()(t,cd,mp);
    }
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask(Task& t, const D& cd, ManagePtr& mp)
    {
    return cloneDoTask_Case1ConstMP<Ret>(t,cd,mp,0);
    }

/////////////////////

//(4) Least preferred version which always compiles
template <typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask_Impl(long l, Task& t, D1& d1, const D2& d2)
    {
    Error("2 parameter doTask not defined for specified task or data type");
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
    -> std::conditional_t<std::is_void<decltype(doTask(t,d1,d2,mp))>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,d1,d2,mp));
    return FixRet<Ret,ActualRet,Task&,D1&,const D2&,ManagePtr&>()(t,d1,d2,mp);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask(Task& t, D1& d1, const D2& d2, ManagePtr& mp)
    {
    //First try calling function labeled (1) above
    return callDoTask_Impl<Ret,Task,D1,D2>(0,t,d1,d2,mp);
    }

/////////////////////

//(2-fail)
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_Case2ConstNoMP(Task& t, const D1& cd1, const D2& d2, ManagePtr& mp,long) 
    {
    if(!mp.parg1().unique()) mp.parg1() = mp.parg1()->clone();
    auto* pd1 = static_cast<ITDataType<std::remove_const_t<D1>>*>(mp.parg1().get());
    return callDoTask<Ret>(t,pd1->d,d2,mp);
    }
//(2-success)
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_Case2ConstNoMP(Task& t, const D1& cd1, const D2& d2, ManagePtr& mp,int) 
    -> std::conditional_t<std::is_same<decltype(doTask(t,cd1,d2)),void>::value,Ret,Ret>
    {
    return FixRet<Ret,decltype(doTask(t,cd1,d2)),Task&,const D1&,const D2&>()(t,cd1,d2);
    }
//(1-fail)
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_Case1ConstMP(Task& t, const D1& cd1, const D2& d2, ManagePtr& mp,long)
    {
    return cloneDoTask_Case2ConstNoMP<Ret,Task,D1,D2>(t,cd1,d2,mp,0);
    }
//(1-success)
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_Case1ConstMP(Task& t, const D1& cd1, const D2& d2, ManagePtr& mp,int) 
    -> std::conditional_t<std::is_same<decltype(doTask(t,cd1,d2,mp)),void>::value,Ret,Ret>
    {
    return FixRet<Ret,decltype(doTask(t,cd1,d2,mp)),Task&,const D1&,const D2&,ManagePtr&>()(t,cd1,d2,mp);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask(Task& t, const D1& d1, const D2& d2, ManagePtr& mp)
    {
    return cloneDoTask_Case1ConstMP<Ret,Task,D1,D2>(t,d1,d2,mp,0);
    }

void inline
check(const CPData& p)
    {
    if(!p) Error("doTask called on unallocated store pointer");
    }

} //namespace detail

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
doTask(Task&& t,
       const ITData& arg)
    {
    RegisterTask<OneArg,Task,ReturnType> r(std::forward<Task>(t));
    arg.plugInto(r);
    return std::move(r.getReturn());
    }

template<typename Task>
Task
doTask(Task&& t,
       const ITData& arg)
    {
    RegisterTask<OneArg,Task,Void> r(std::forward<Task>(t));
    arg.plugInto(r);
    return std::move(r.getTask());
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task&& t,
       const CPData& arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    RegisterTask<OneArg,Task,ReturnType> r(std::forward<Task>(t));
    arg->plugInto(r);
    return std::move(r.getReturn());
    }
template<typename Task>
Task
doTask(Task&& t,
       const CPData& arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    RegisterTask<OneArg,Task,Void> r(std::forward<Task>(t));
    arg->plugInto(r);
    return std::move(r.getTask());
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task&& t,
       PData& arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    RegisterTask<OneArg,Task,ReturnType> r(std::forward<Task>(t),&arg);
    arg->plugInto(r);
    return std::move(r.getReturn());
    }
template<typename Task>
Task
doTask(Task&& t,
       PData& arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    RegisterTask<OneArg,Task,Void> r(std::forward<Task>(t),&arg);
    arg->plugInto(r);
    return std::move(r.getTask());
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task&& t,
       const CPData& arg1,
       const CPData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    RegisterTask<TwoArgs,Task,ReturnType> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getReturn());
    }
template<typename Task>
Task
doTask(Task&& t,
       const CPData& arg1,
       const CPData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    RegisterTask<TwoArgs,Task,Void> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getTask());
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task&& t,
       PData& arg1,
       const CPData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    RegisterTask<TwoArgs,Task,ReturnType> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getReturn());
    }
template<typename Task>
Task
doTask(Task&& t,
       PData& arg1,
       const CPData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    RegisterTask<TwoArgs,Task,Void> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getTask());
    }

template<typename ReturnType, typename Task>
ReturnType
doTask(Task&& t,
       PData& arg1,
       const ITData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
#endif
    RegisterTask<TwoArgs,Task,ReturnType> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getReturn());
    }
template<typename Task>
Task
doTask(Task&& t,
       PData& arg1,
       const ITData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
#endif
    RegisterTask<TwoArgs,Task,Void> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getTask());
    }

} //namespace itensor

#undef LPAREN
#undef RPAREN
#undef REGISTER_TYPES

#endif

