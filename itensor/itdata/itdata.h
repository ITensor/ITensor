//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_H
#define __ITENSOR_ITDATA_H
#include "../global.h"
#include "../detail/algs.h"
#include "../detail/call_rewrite.h"

//
// To register a new ITData subtype:
//
// (1) Add a new line to the REGISTER macro below, following the same format
//     and no trailing \ on the last line.
// (2) Forward-declare the subtype just after "namespace itensor" below.
//

////////////////////////////////////
// (1) Add a new line here to register a new ITData subtype:
#define REGISTER(X,Y,Z)\
    X(Y  ITDense<Real>      &d) Z\
    X(Y  ITDense<Complex>   &d) Z\
    X(Y  ITCombiner         &d) Z\
    X(Y  ITDiag<Real>       &d) Z\
    X(Y  ITDiag<Complex>    &d) Z\
    X(Y  IQTData<Real>      &d) Z
///////////////////////////////////

namespace itensor {

///////////////////////////////////
// (2) Forward declarations for all ITData subtypes:
template <typename>
class ITDense;

template <typename>
class ITDiag;

class ITCombiner;

template <typename>
class IQTData;
///////////////////////////////////


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

    REGISTER(void virtual applyTo,,=0;)
    REGISTER(void virtual applyTo,const,=0;)

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

template <class Derived>
struct RegisterData : ITData
    {
    RegisterData() { }
    virtual ~RegisterData() { }

    private:
    
    NewData
    clone() const final 
        { 
        auto* pdt = static_cast<const Derived*>(this);
        return std::make_unique<Derived>(*pdt);
        }

    void
    plugInto(FuncBase& f) const final
        {
        auto& cdt = *(static_cast<const Derived*>(this));
        f.applyTo(cdt);
        }
        
    void
    plugInto(FuncBase& f) final
        {
        auto& dt = *(static_cast<Derived*>(this));
        f.applyTo(dt);
        }
    };

//////////////////
//////////////////

struct NoReturn { };

template <typename Derived, typename Return = NoReturn>
struct RegisterFunc : FuncBase
    {
    using return_type = Return;

    RegisterFunc();
    RegisterFunc(RegisterFunc&& other);
    virtual ~RegisterFunc() { }

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

    operator Return() { return ret_; }

    private:

    void
    updateArg1();

    template<typename T>
    void
    applyToImpl(const T& d);

    template<typename T>
    void
    applyToImpl(T& d);

    enum Action
        {
        None,
        AssignNewData,
        AssignPointerRtoL
        };

    PData* arg1_ = nullptr;
    const PData* arg2_ = nullptr;
    Derived& dt_;
    Action action_ = None;
    PData nd_;
    Return ret_;

    public:

    void
    setup(PData* arg1);
    void
    setup(PData* arg1, const PData* arg2);

    REGISTER(void applyTo,,final { applyToImpl(d); })
    REGISTER(void applyTo,const,final { applyToImpl(d); })

    };

template <typename Callable, typename T1, typename Return>
struct CallWrap : FuncBase
    {
    CallWrap(Callable& c, T1& arg1) : c_(c), arg1_(arg1), parg1_(nullptr) { }
    CallWrap(Callable& c, T1& arg1, PData& parg1) : c_(c), arg1_(arg1), parg1_(&parg1) { }

    REGISTER(void applyTo,,final { applyToImpl(d); })
    REGISTER(void applyTo,const,final { applyToImpl(d); })

    Return
    getReturn() { return ret_; }

    private:

    template<typename T2>
    void
    applyToImpl(const T2& d2);

    Callable& c_;
    T1& arg1_;
    PData* parg1_;
    Return ret_;
    };


//
// Implementations
//

namespace detail {

template<typename Ret, typename Func, typename T>
auto
clone_modify_impl(Func& f, T& a, PData& pdat,int) //-> decltype(f(static_cast<const T&>(a)))
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T&) to get substitution failure (SFINAE) if no such call exists.
    -> std::conditional_t<std::is_same<decltype(f(static_cast<const T&>(a))),void>::value,Ret,Ret>
    {
    const T& ca = a;
    return detail::call<Ret>(f,ca);
    }

template<typename Ret, typename Func, typename T>
Ret
clone_modify_impl(Func& f, T& a, PData& pdat,long)
    {
    //if(Global::debug3()) println("Calling solo (1 param)");
    //println("--> Calling solo (1 param)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pa = static_cast<T*>(pdat.get());
    return detail::call<Ret>(f,*pa);
    }

template<typename Ret, typename Func, typename T>
Ret
clone_modify(Func& f, T& a, PData& pdat)
    {
    return clone_modify_impl<Ret>(f,a,pdat,0);
    }

/////////////////////

template<typename Ret, typename Func, typename T1, typename T2>
auto
clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,int) 
    //Using std::conditional here because we want the return type to be Ret regardless
    //but need to call possibly void f(const T1&,const T2&) to get substitution failure (SFINAE) 
    //if no such call exists.
    -> std::conditional_t<std::is_same<decltype(f(static_cast<const T1&>(a1),a2)),void>::value,Ret,Ret>
    {
    //if(Global::debug3()) println("Not calling solo (2 params)");
    const T1& ca1 = a1;
    return detail::call<Ret>(f,ca1,a2);
    }

template<typename Ret, typename Func, typename T1, typename T2>
Ret
clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,long)
    {
    //if(Global::debug3()) println("Calling solo (2 params)");
    //println("--> Calling solo (2 params)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pa1 = static_cast<T1*>(pdat.get());
    return detail::call<Ret>(f,*pa1,a2);
    }

template<typename Ret, typename Func, typename T1, typename T2>
Ret
clone_modify(Func& f, T1& a1, const T2& a2, PData& pdat)
    {
    return clone_modify_impl<Ret>(f,a1,a2,pdat,0);
    }

}; //namespace detail

template <typename Callable, typename T1, typename Return>
template<typename T2>
void CallWrap<Callable,T1,Return>::
applyToImpl(const T2& d2) 
    { 
    if(parg1_) ret_ = detail::clone_modify<Return>(c_,arg1_,d2,*parg1_); 
    else       ret_ = detail::call<Return>(c_,arg1_,d2); 
    }

template <typename Derived, typename Return>
RegisterFunc<Derived,Return>::
RegisterFunc() 
    :
    arg2_(nullptr),
    dt_(*static_cast<Derived*>(this))
    { }

template <typename Derived, typename Return>
RegisterFunc<Derived,Return>::
RegisterFunc(RegisterFunc&& other)
    :
    arg1_(other.arg1_),
    arg2_(other.arg2_),
    dt_(*static_cast<Derived*>(this)),
    action_(other.action_),
    nd_(std::move(other.nd_)),
    ret_(std::move(other.ret_))
    { 
    other.arg1_ = nullptr;
    other.arg2_ = nullptr;
    other.action_ = None;
    }

template <typename Derived, typename Return>
void RegisterFunc<Derived,Return>::
setup(PData* arg1)
    { 
    arg1_ = arg1;
    }

template <typename Derived, typename Return>
void RegisterFunc<Derived,Return>::
setup(PData* arg1, const PData* arg2)
    { 
    arg1_ = arg1; 
    arg2_ = arg2;
    }

template <typename Derived, typename Return>
template <typename ITDataType, typename... Args>
ITDataType* RegisterFunc<Derived,Return>::
makeNewData(Args&&... args)
    {
    if(!arg1_) Error("Can't call setNewData with const-only access to first arg");
    action_ = AssignNewData;
    auto newdat = std::make_shared<ITDataType>(std::forward<Args>(args)...);
    auto* ret = newdat.get();
    nd_ = std::move(newdat);
    return ret;
    }

template <typename Derived, typename Return>
template <typename ITDataType>
void RegisterFunc<Derived,Return>::
setNewData(std::shared_ptr<ITDataType>&& nd)
    {
    nd_ = std::move(nd);
    action_ = AssignNewData;
    }

template <typename Derived, typename Return>
void RegisterFunc<Derived,Return>::
assignPointerRtoL() 
    { 
    if(!arg2_) Error("No second pointer provided for action AssignPointerRtoL");
    action_ = AssignPointerRtoL; 
    }

template <typename Derived, typename Return>
void RegisterFunc<Derived,Return>::
updateArg1()
    {
    if(action_ == AssignNewData)
        {
        *arg1_ = std::move(nd_);
        }
    else if(action_ == AssignPointerRtoL)
        {
        *arg1_ = *arg2_;
        }
    }

template <typename Derived, typename Return>
template<typename T>
T& RegisterFunc<Derived,Return>::
modifyData(const T& d)
    {
    if(!arg1_) Error("Can't modify const data");
    if(!(arg1_->unique())) *arg1_ = (*arg1_)->clone();
    auto* pa1 = static_cast<T*>(arg1_->get());
    return *pa1;
    }

template <typename Derived, typename Return>
template<typename T>
void RegisterFunc<Derived,Return>::
applyToImpl(const T& d)
    {
    //println("In applyToImpl #1");
    if(arg2_)
        {
        CallWrap<Derived,const T,Return> w(dt_,d);
        (*arg2_)->plugInto(w);
        ret_ = std::move(w.getReturn());
        }
    else
        {
        ret_ = detail::call<Return>(dt_,d);
        }
    updateArg1();
    }

template <typename Derived, typename Return>
template<typename T>
void RegisterFunc<Derived,Return>::
applyToImpl(T& d)
    {
    //println("In applyToImpl #2");
    assert(arg1_);
    if(arg2_)
        {
        CallWrap<Derived,T,Return> w(dt_,d,*arg1_);
        (*arg2_)->plugInto(w);
        ret_ = std::move(w.getReturn());
        }
    else
        {
        ret_ = detail::clone_modify<Return>(dt_,d,*arg1_);
        }
    updateArg1();
    }


//////
////// applyFunc methods
//////

template <typename F>
using returnTypeOf = typename std::conditional<
                               std::is_same<typename F::return_type,NoReturn>::value,
                                   F,
                                   typename F::return_type
                               >::type;

template<typename F>
returnTypeOf<F>
applyFunc(const ITData& arg,
          F f = F())
    {
    arg.plugInto(f);
    return f;
    }

template<typename F>
returnTypeOf<F>
applyFunc(const CPData& arg,
          F f = F())
    {
    arg->plugInto(f);
    return f;
    }

template<typename F>
returnTypeOf<F>
applyFunc(PData& arg,
          F f = F())
    {
    f.setup(&arg);
    arg->plugInto(f);
    return f;
    }

template<typename F>
returnTypeOf<F>
applyFunc(const PData& arg1,
          const PData& arg2,
          F f = F())
    {
    f.setup(&arg1,&arg2);
    arg2->plugInto(f);
    return f;
    }

template<typename F>
returnTypeOf<F>
applyFunc(PData& arg1,
          const PData& arg2,
          F f = F())
    {
    f.setup(&arg1,&arg2);
    arg1->plugInto(f);
    return f;
    }

}; //namespace itensor

#endif

