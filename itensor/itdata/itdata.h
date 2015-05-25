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
using NewData = std::unique_ptr<ITData>;

template<typename DataType, typename... Args>
std::unique_ptr<DataType>
make_newdata(Args&&... args)
    {
    return std::unique_ptr<DataType>(new DataType(std::forward<Args>(args)...));
    }


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

    NewData virtual
    clone() const = 0;

    void virtual
    plugInto(FuncBase& f) const = 0;

    void virtual
    plugInto(FuncBase& f) = 0;
    };

template <class DType>
struct ITDataType : ITData
    {
    public:
    DType d;

    ITDataType() { }

    template<typename... VArgs>
    ITDataType(VArgs&&... vargs)
        : d(std::forward<VArgs>(vargs)...)
        { }

    virtual ~ITDataType() { }

    private:
    
    NewData
    clone() const final 
        { 
        return std::make_unique<ITDataType<DType>>(d);
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


    template <typename T>
    T&
    modifyData(const T& d);

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

template <typename Task, typename Return = Void>
class RegisterTask : public FuncBase
    {
    private:
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
        ret_(std::move(o.ret_)),
        mp_(std::move(o.mp_))
        { }

    virtual ~RegisterTask() { }

    operator Return() { return std::move(ret_); }

    private:

    template<typename D>
    void
    applyToImpl(const D& d);

    template<typename D>
    void
    applyToImpl(D& d);


    public:

    REGISTER_TYPES(void applyTo LPAREN, &d RPAREN final { applyToImpl(d); } )
    REGISTER_TYPES(void applyTo LPAREN, const&d RPAREN final { applyToImpl(d); } )

    };

//template <typename Callable, typename T1, typename Return>
//struct CallWrap : FuncBase
//    {
//    CallWrap(Callable& c, T1& arg1) : c_(c), arg1_(arg1), parg1_(nullptr) { }
//    CallWrap(Callable& c, T1& arg1, PData& parg1) : c_(c), arg1_(arg1), parg1_(&parg1) { }
//
//    REGISTER(void applyTo,,final { applyToImpl(d); })
//    REGISTER(void applyTo,const,final { applyToImpl(d); })
//
//    Return
//    getReturn() { return ret_; }
//
//    private:
//
//    template<typename D2>
//    void
//    applyToImpl(const D2& d2)
//        { 
//        if(parg1_) ret_ = detail::cloneDoTask<Return>(c_,arg1_,d2,*parg1_); 
//        else       ret_ = detail::call<Return>(c_,arg1_,d2); 
//        }
//    
//    Callable& c_;
//    T1& arg1_;
//    PData* parg1_;
//    Return ret_;
//    };



//
// Implementations
//

namespace detail {

//If ActualRet!=void this gets called
template<typename Ret, typename Task, typename D, typename ActualRet>
struct FixRet
    {
    Ret
    operator()(Task& t, D& d)
        {
        static_assert(!std::is_same<Ret,Void>::value,
                      "No return type specified in call to doTask");
        static_assert(std::is_same<Ret,ActualRet>::value,
                      "Mismatched return type specified in call to doTask");
        return doTask(t,d);
        }
    };
//If ActualRet==void this gets called
template<typename Ret, typename Task, typename D>
struct FixRet<Ret,Task,D,void>
    {
    Ret
    operator()(Task& t, D& d)
        {
        doTask(t,d);
        return Ret{};
        }
    };
//Less preferred version which always compiles
template <typename Ret, typename Task, typename D>
Ret
CallDoTask_Impl(Task& t, D& d, long) 
    {
    throw std::runtime_error("doTask not defined for specified task or data type");
    return Ret{};
    }
//This is the preferred version since 0->int requires no cast
//Template substitution will fail (not a compile-time error)
//if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret, typename Task, typename D>
auto 
CallDoTask_Impl(Task& t, D& d, int) 
    -> std::conditional_t<std::is_same<decltype(doTask(t,d)),void>::value,Ret,Ret>
    {
    using ActualRet = decltype(doTask(t,d));
    return FixRet<Ret,Task,D,ActualRet>()(t,d);
    }
//CallDoTask attempts the following: return doTask(t,d);
//- If doTask(t,d) not defined, throws an exception
//  if called (this converts what would otherwise be 
//  a compile-time error into a run-time error)
//- If doTask(t,d) defined but returns void, 
//  returns Ret{} instead
template<typename Ret, typename Task, typename D>
Ret
CallDoTask(Task& t, D& d)
    {
    return CallDoTask_Impl<Ret,Task,D>(t,d,0);
    }

//Less preferred version calls doTask(Task,D&)
//and makes sure pdat is unique first ("copy on write")
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Impl(Task& t, D& d, PData& pdat,long)
    {
    //if(Global::debug3()) println("Calling solo (1 param)");
    //println("--> Calling solo (1 param)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pd = static_cast<ITDataType<D>*>(pdat.get());
    return CallDoTask<Ret>(t,pd->d);
    }
//Preferred version since 0->int requires no conversion
//Attempts to call doTask(Task,const D&) if defined
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Impl(Task& t, D& d, PData& pdat,int)
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T&) to get substitution failure (SFINAE) if no such call exists.
    -> std::conditional_t<std::is_same<decltype(doTask(t,static_cast<const D&>(d))),void>::value,Ret,Ret>
    {
    //println("--> Not calling solo (1 param)");
    const auto& cd = d;
    return CallDoTask<Ret>(t,cd);
    }
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask(Task& t, D& d, PData& pdat)
    {
    return cloneDoTask_Impl<Ret>(t,d,pdat,0);
    }

/////////////////////

//template<typename Ret, typename Func, typename T1, typename T2>
//auto
//clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,int) 
//    //Using std::conditional here because we want the return type to be Ret regardless
//    //but need to call possibly void f(const T1&,const T2&) to get substitution failure (SFINAE) 
//    //if no such call exists.
//    -> std::conditional_t<std::is_same<decltype(f(static_cast<const T1&>(a1),a2)),void>::value,Ret,Ret>
//    {
//    //if(Global::debug3()) println("Not calling solo (2 params)");
//    const T1& ca1 = a1;
//    return detail::call<Ret>(f,ca1,a2);
//    }
//
//template<typename Ret, typename Func, typename T1, typename T2>
//Ret
//clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,long)
//    {
//    //if(Global::debug3()) println("Calling solo (2 params)");
//    //println("--> Calling solo (2 params)");
//    if(!pdat.unique()) pdat = pdat->clone();
//    auto* pa1 = static_cast<T1*>(pdat.get());
//    return detail::call<Ret>(f,*pa1,a2);
//    }
//
//template<typename Ret, typename Func, typename T1, typename T2>
//Ret
//clone_modify(Func& f, T1& a1, const T2& a2, PData& pdat)
//    {
//    return clone_modify_impl<Ret>(f,a1,a2,pdat,0);
//    }

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
    return ret;
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
T& ManagePtr::
modifyData(const T& d)
    {
    if(!parg1_) Error("Can't modify const data");
    if(!(parg1_->unique())) 
        {
        *parg1_ = (*parg1_)->clone();
        }
    auto* pa1 = static_cast<T*>(parg1_->get());
    return *pa1;
    }

template<typename Task, typename Ret>
template<typename D>
void RegisterTask<Task,Ret>::
applyToImpl(const D& d)
    {
    //println("In applyToImpl #1");
    if(mp_.hasArg2())
        {
        throw std::runtime_error("Two-arg doTask not yet implemented");
        //CallWrap<Task,const D,Return> w(t_,d);
        //mp_.arg2()->plugInto(w);
        //ret_ = std::move(w.getReturn());
        }
    else
        {
        ret_ = detail::CallDoTask<Ret>(t_,d);
        }
    }

template<typename Task, typename Ret>
template<typename D>
void RegisterTask<Task,Ret>::
applyToImpl(D& d)
    {
    //println("In applyToImpl #2");
    assert(mp_.hasPArg1());
    if(mp_.hasArg2())
        {
        throw std::runtime_error("Two-arg doTask not yet implemented");
        //CallWrap<Task,D,Return> w(t_,d,mp.parg1());
        //mp_.arg2()->plugInto(w);
        //ret_ = std::move(w.getReturn());
        }
    else
        {
        ret_ = detail::cloneDoTask<Ret>(t_,d,mp_.parg1());
        }
    }


//////
////// doTask methods
//////


template<typename ReturnType, typename Task>
ReturnType
doTask(Task t,
       const ITData& arg)
    {
    RegisterTask<Task,ReturnType> r(std::move(t));
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
    RegisterTask<Task,ReturnType> r(std::move(t));
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
    RegisterTask<Task,ReturnType> r(std::move(t),&arg);
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

//template<typename Ret, typename Task>
//Ret
//doTaskReturn(Task t,
//             const PData& arg1,
//             const PData& arg2)
//    {
//    RegisterTask<Task,Ret> r(std::move(t),&arg1,&arg2);
//    arg1->plugInto(r);
//    return r;
//    }
//template<typename Task>
//Task
//doTask(Task&& t,
//       const PData& arg1,
//       const PData& arg2)
//    {
//    doTaskReturn<Void,Task>(t,arg1,arg2);
//    return t;
//    }
//
//template<typename Ret, typename Task>
//Ret
//doTaskReturn(Task t,
//             PData& arg1,
//             const PData& arg2)
//    {
//    RegisterTask<Task,Ret> r(std::move(t),&arg1,&arg2);
//    arg1->plugInto(r);
//    return r;
//    }
//template<typename Task>
//Task
//doTask(Task&& t,
//       PData& arg1,
//       const PData& arg2)
//    {
//    doTaskReturn<Void,Task>(t,arg1,arg2);
//    return t;
//    }
//
//template<typename Ret, typename Task>
//Ret
//doTaskReturn(Task t,
//             PData& arg1,
//             const ITData& arg2)
//    {
//    RegisterTask<Task,Ret> r(std::move(t),&arg1,&arg2);
//    arg1->plugInto(r);
//    return r;
//    }
//template<typename Task>
//Task
//doTask(Task&& t,
//       PData& arg1,
//       const ITData& arg2)
//    {
//    doTaskReturn<Void,Task>(t,arg1,arg2);
//    return t;
//    }

} //namespace itensor

#undef LPAREN
#undef RPAREN
#undef REGISTER_TYPES

#endif

