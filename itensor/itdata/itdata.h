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
        return f.applyTo(cdt);
        }
        
    void
    plugInto(FuncBase& f) final
        {
        auto& dt = *(static_cast<Derived*>(this));
        return f.applyTo(dt);
        }
    };

//////////////////
//////////////////


template <typename Derived>
struct RegisterFunc : FuncBase
    {
    RegisterFunc();
    RegisterFunc(RegisterFunc&& other);
    virtual ~RegisterFunc() { }

    template <typename T>
    T&
    modifyData(const T& d);

    template <typename ITDataType, typename... Args>
    ITDataType*
    setNewData(Args&&... args);

    bool
    newDataIsSet() { return action_ == AssignNewData; }

    void
    assignPointerRtoL();

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

    public:

    void
    setup(PData* arg1);
    void
    setup(PData* arg1, const PData* arg2);

    REGISTER(void applyTo,,final { applyToImpl(d); })
    REGISTER(void applyTo,const,final { applyToImpl(d); })

    };

template <typename Callable, typename T1>
struct CallWrap : FuncBase
    {
    CallWrap(Callable& c, T1& arg1) : c_(c), arg1_(arg1), parg1_(nullptr) { }
    CallWrap(Callable& c, T1& arg1, PData& parg1) : c_(c), arg1_(arg1), parg1_(&parg1) { }

    REGISTER(void applyTo,,final { applyToImpl(d); })
    REGISTER(void applyTo,const,final { applyToImpl(d); })

    private:

    template<typename T2>
    void
    applyToImpl(const T2& d2);

    Callable& c_;
    T1& arg1_;
    PData* parg1_;
    };


//
// Implementations
//

namespace detail {

template<typename Func, typename T>
auto
clone_modify_impl(Func& f, T& a, PData& pdat,int) //-> decltype(f(static_cast<const T&>(a)))
    //Using std::conditional here because we want the return type to be a placeholder int,
    //but need to call f(const T&) to get substitution failure (SFINAE) if no such call exists.
    -> std::conditional_t<std::is_same<decltype(f(static_cast<const T&>(a))),void>::value,int,int>
    {
    const T& ca = a;
    detail::call(f,ca);
    return 0;
    }

template<typename Func, typename T>
void
clone_modify_impl(Func& f, T& a, PData& pdat,long)
    {
    //if(Global::debug3()) println("Calling solo (1 param)");
    //println("--> Calling solo (1 param)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pa = static_cast<T*>(pdat.get());
    detail::call(f,*pa);
    }

template<typename Func, typename T>
void
clone_modify(Func& f, T& a, PData& pdat)
    {
    clone_modify_impl(f,a,pdat,0);
    }

/////////////////////

template<typename Func, typename T1, typename T2>
auto
clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,int) 
    //Using std::conditional here because we want the return type to be a placeholder int,
    //but need to call f(const T1&,const T2&) to get substitution failure (SFINAE) 
    //if no such call exists.
    -> std::conditional_t<std::is_same<decltype(f(static_cast<const T1&>(a1),a2)),void>::value,int,int>
    {
    //if(Global::debug3()) println("Not calling solo (2 params)");
    const T1& ca1 = a1;
    detail::call(f,ca1,a2);
    return 0;
    }

template<typename Func, typename T1, typename T2>
void
clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,long)
    {
    //if(Global::debug3()) println("Calling solo (2 params)");
    //println("--> Calling solo (2 params)");
    if(!pdat.unique()) pdat = pdat->clone();
    auto* pa1 = static_cast<T1*>(pdat.get());
    detail::call(f,*pa1,a2);
    }

template<typename Func, typename T1, typename T2>
void
clone_modify(Func& f, T1& a1, const T2& a2, PData& pdat)
    {
    clone_modify_impl(f,a1,a2,pdat,0);
    }

}; //namespace detail

template <typename Callable, typename T1>
template<typename T2>
void CallWrap<Callable,T1>::
applyToImpl(const T2& d2) 
    { 
    if(parg1_) detail::clone_modify(c_,arg1_,d2,*parg1_); 
    else       detail::call(c_,arg1_,d2); 
    }

template <typename Derived>
RegisterFunc<Derived>::
RegisterFunc() 
    :
    arg2_(nullptr),
    dt_(*static_cast<Derived*>(this))
    { }

template <typename Derived>
RegisterFunc<Derived>::
RegisterFunc(RegisterFunc&& other)
    :
    arg1_(other.arg1_),
    arg2_(other.arg2_),
    dt_(*static_cast<Derived*>(this)),
    action_(other.action_),
    nd_(std::move(other.nd_))
    { 
    other.arg1_ = nullptr;
    other.arg2_ = nullptr;
    other.action_ = None;
    }

template <typename Derived>
void RegisterFunc<Derived>::
setup(PData* arg1)
    { 
    arg1_ = arg1;
    }

template <typename Derived>
void RegisterFunc<Derived>::
setup(PData* arg1, const PData* arg2)
    { 
    arg1_ = arg1; 
    arg2_ = arg2;
    }

template <typename Derived>
template <typename ITDataType, typename... Args>
ITDataType* RegisterFunc<Derived>::
setNewData(Args&&... args)
    {
    if(!arg1_) Error("Can't call setNewData with const-only access to first arg");
    action_ = AssignNewData;
    auto newdat = std::make_shared<ITDataType>(std::forward<Args>(args)...);
    auto* ret = newdat.get();
    nd_ = std::move(newdat);
    return ret;
    }

template <typename Derived>
void RegisterFunc<Derived>::
assignPointerRtoL() 
    { 
    if(!arg2_) Error("No second pointer provided for action AssignPointerRtoL");
    action_ = AssignPointerRtoL; 
    }

template <typename Derived>
void RegisterFunc<Derived>::
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

template <typename Derived>
template<typename T>
T& RegisterFunc<Derived>::
modifyData(const T& d)
    {
    if(!arg1_) Error("Can't modify const data");
    if(!(arg1_->unique())) *arg1_ = (*arg1_)->clone();
    auto* pa1 = static_cast<T*>(arg1_->get());
    return *pa1;
    }

template <typename Derived>
template<typename T>
void RegisterFunc<Derived>::
applyToImpl(const T& d)
    {
    //println("In applyToImpl #1");
    if(arg2_)
        {
        CallWrap<Derived,const T> w(dt_,d);
        (*arg2_)->plugInto(w);
        }
    else
        {
        detail::call(dt_,d);
        }
    updateArg1();
    }

template <typename Derived>
template<typename T>
void RegisterFunc<Derived>::
applyToImpl(T& d)
    {
    //println("In applyToImpl #2");
    assert(arg1_);
    if(arg2_)
        {
        CallWrap<Derived,T> w(dt_,d,*arg1_);
        (*arg2_)->plugInto(w);
        }
    else
        {
        detail::clone_modify(dt_,d,*arg1_);
        }
    updateArg1();
    }


//////
////// applyFunc methods
//////

template<typename F>
F
applyFunc(const ITData& arg,
          F f = F())
    {
    arg.plugInto(f);
    return f;
    }

template<typename F>
F
applyFunc(const CPData& arg,
          F f = F())
    {
    arg->plugInto(f);
    return f;
    }

template<typename F>
F
applyFunc(PData& arg,
          F f = F())
    {
    f.setup(&arg);
    arg->plugInto(f);
    return f;
    }

template<typename F>
F
applyFunc(const PData& arg1,
          const PData& arg2,
          F f = F())
    {
    f.setup(&arg1,&arg2);
    arg2->plugInto(f);
    return f;
    }

template<typename F>
F 
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

