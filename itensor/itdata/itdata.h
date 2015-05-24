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
    DType d_;
    public:

    ITDataType() { }

    template<typename... VArgs>
    ITDataType(VArgs&&... vargs)
        : d_(std::forward<VArgs>(vargs)...)
        { }

    virtual ~ITDataType() { }

    private:
    
    NewData
    clone() const final 
        { 
        return std::make_unique<ITDataType<DType>>(d_);
        }

    void
    plugInto(FuncBase& f) const final
        {
        f.applyTo(d_);
        }
        
    void
    plugInto(FuncBase& f) final
        {
        f.applyTo(d_);
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

    template <typename ITDataType, typename... Args>
    ITDataType*
    makeNewData(Args&&... args);

    template <typename ITDataType>
    void
    setNewData(std::shared_ptr<ITDataType>&& nd);

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

template<typename T1, typename T2>
void
doTask(const T1& t1, T2& t2)
    {
    throw std::runtime_error("doTask not defined for specified task or data type");
    }

namespace detail {

//template <typename Ret, class Task, typename D>
//auto 
//callDoTask_Impl(Task&& t, D&& d, int) 
//    -> decltype(doTask(std::forward<Task>(t),std::forward<D>(d)))
//    {
//    return doTask(std::forward<Task>(t),std::forward<D>(d));
//    }
//template <typename Ret, class Task, typename D>
//Ret
//callDoTask_Impl(Task&& t, D&& d, long) 
//    {
//    throw std::runtime_error("doTask not defined for specified task or data type");
//    return Ret{};
//    }

//template <typename Ret, class Task, typename D>
//auto 
//callDoTask_Impl(Task&& t, D&& d, int) 
//    //-> std::conditional_t<std::is_same<decltype(doTask(std::forward<Task>(t),std::forward<D>(d))),void>::value,Ret,Ret>
//    -> std::conditional_t<std::is_same<decltype(doTask(t,std::forward<D>(d))),void>::value,Ret,Ret>
//    {
//    doTask(std::forward<Task>(t),std::forward<D>(d));
//    return Ret{};
//    }
//template <typename Ret, class Task, typename D>
//Ret
//callDoTask_Impl(Task&& t, D&& d, long) 
//    {
//    throw std::runtime_error("doTask not defined for specified task or data type");
//    return Ret{};
//    }

//template<typename Ret, typename Task, typename D, typename ActualRet>
//struct FixRet
//    {
//    Ret
//    operator()(Task&& t, D&& d) const
//        {
//        return callDoTask_Impl<Ret,Task,D>(std::forward<Task>(t),std::forward<D>(d),0);
//        }
//    };
//template<typename Ret, typename Task, typename D>
//struct FixRet<Ret,Task,D,void> //specialization for case ActualRet=void
//    {
//    Ret
//    operator()(Task&& t, D&& d) const
//        {
//        callDoTask_Impl<void,Task,D>(std::forward<Task>(t),std::forward<D>(d),0);
//        return Ret{};
//        }
//    };
template <typename Ret, class Task, typename D>
Ret
callDoTask(Task&& t, D&& d)
    {
    //using ActualRet = std::result_of_t<decltype(doTask)&(Task,D)>;
    //return FixRet<Ret,Task,D,ActualRet>()(std::forward<Task>(t),std::forward<D>(d),0);
    //return callDoTask_Impl<Ret,Task,D>(std::forward<Task>(t),std::forward<D>(d),0);
    doTask(std::forward<Task>(t),std::forward<D>(d));
    return Ret{};
    }

template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Impl(Task& t, D& d, PData& pdat,int)
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T&) to get substitution failure (SFINAE) if no such call exists.
    -> std::conditional_t<std::is_same<decltype(doTask(t,static_cast<const D&>(d))),void>::value,Ret,Ret>
    {
    const auto& cd = d;
    return callDoTask<Ret>(t,cd);
    }
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Impl(Task& t, D& d, PData& pdat,long)
    {
    //if(Global::debug3()) println("Calling solo (1 param)");
    println("--> Calling solo (1 param)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pd = static_cast<D*>(pdat.get());
    return callDoTask<Ret>(t,*pd);
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

template <typename ITDataType, typename... VArgs>
ITDataType* ManagePtr::
makeNewData(VArgs&&... vargs)
    {
    if(!parg1_) Error("Can't call makeNewData with const-only access to first arg");
    action_ = AssignNewData;
    auto newdat = std::make_shared<ITDataType>(std::forward<VArgs>(vargs)...);
    auto* ret = newdat.get();
    nd_ = std::move(newdat);
    return ret;
    }

template <typename ITDataType>
void ManagePtr::
setNewData(std::shared_ptr<ITDataType>&& nd)
    {
    nd_ = std::move(nd);
    action_ = AssignNewData;
    }

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
        ret_ = detail::callDoTask<Ret>(t_,d);
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


template<typename Ret, typename Task>
Ret
doTaskReturn(Task t,
             const ITData& arg)
    {
    RegisterTask<std::remove_reference_t<Task>,Ret> r(std::move(t));
    arg.plugInto(r);
    return r;
    }
template<typename Task>
Task
doTask(Task&& t,
       const ITData& arg)
    {
    doTaskReturn<Void,Task>(std::forward<Task>(t),arg);
    return t;
    }

//template<typename Ret, typename Task>
//Ret
//doTaskReturn(Task t,
//             const CPData& arg)
//    {
//    RegisterTask<Task,Ret> r(std::move(t));
//    arg->plugInto(r);
//    return r;
//    }
//template<typename Task>
//Task
//doTask(Task&& t,
//       const CPData& arg)
//    {
//    doTaskReturn<Void,Task>(t,arg);
//    return t;
//    }
//
//template<typename Ret, typename Task>
//Ret
//doTaskReturn(Task t,
//             PData& arg)
//    {
//    RegisterTask<Task,Ret> r(std::move(t),&arg);
//    arg->plugInto(r);
//    return r;
//    }
//template<typename Task>
//Task
//doTask(Task&& t,
//       PData& arg)
//    {
//    doTaskReturn<Void,Task>(t,arg);
//    return t;
//    }
//
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

#endif

