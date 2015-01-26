//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_H
#define __ITENSOR_ITDATA_H
#include "../global.h"
#include "../detail/functions.h"

//
// To register a new ITData subtype:
//
// (1) Add a new line to the REGISTER macro below, following the same format
//     no trailing \ on the last line.
// (2) Forward-declare the subtype just after "namespace itensor" below.
//

////////////////////////////////////
// (1) Add a new line here to register a new ITData subtype:
#define REGISTER(X,Y,Z)\
    X(Y  ITDense<Real>      &t) Z;\
    X(Y  ITDense<Complex>   &t) Z;\
    X(Y  ITCombiner         &t) Z;\
    X(Y  ITDiag<Real>       &t) Z;\
    X(Y  ITDiag<Complex>    &t) Z;\
    X(Y  IQTData<Real>      &t) Z;
///////////////////////////////////

//
//Ideas for improvement:
// o Create SFINAE overloads of applyFunc which check if 
//   functions return void, and plug them into wrappers which
//   don't require that they return a ITResult. May conflict
//   with existing definitions of ITData::plugInto though?
//   Not if these wrappers return ITResult() after calling 
//   the wrapped function.
//

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


class ITData;
using PData = std::shared_ptr<ITData>;
using CPData = std::shared_ptr<const ITData>;
using NewData = std::unique_ptr<ITData>;

struct ITResult
    {
    enum Action
        {
        None,
        AssignNewData,
        AssignPointer
        };
    private:
    NewData nd_;
    Action action_ = None;
    public:

    ITResult(Action act = None)
        :
        action_(act)
        { }

    template<typename ITDataType>
    ITResult(std::unique_ptr<ITDataType>&& nd)
        :
        nd_(std::move(nd)),
        action_(None)
        { 
        if(nd_) action_= AssignNewData;
        }

    void
    update(PData& arg) { update_single(arg); }
    void
    update(NewData& arg) { update_single(arg); }

    void
    update(PData& arg1, const PData& arg2)
        {
        if(action_ == AssignNewData)
            {
            arg1 = std::move(nd_);
            }
        else if(action_ == AssignPointer)
            {
            arg1 = arg2;
            }
        }

    private:
    template<typename T>
    void
    update_single(T& arg)
        {
        if(action_ == AssignNewData) arg = std::move(nd_);
#ifdef DEBUG
        else if(action_ == AssignPointer)
            Error("Can't do AssignPointer action on single pointer.");
#endif
        }
    };

template<typename DataType, typename... Args>
std::unique_ptr<DataType>
make_newdata(Args&&... args)
    {
    return std::unique_ptr<DataType>(new DataType(std::forward<Args>(args)...));
    }

template<typename DataType, typename... Args>
ITResult
make_result(Args&&... args)
    {
    return ITResult(std::unique_ptr<DataType>(new DataType(std::forward<Args>(args)...)));
    }

namespace detail {

template<typename Func, typename T>
ITResult
clone_modify(Func& f, T& a, PData& pdat);

template<typename Func, typename T1, typename T2>
ITResult
clone_modify(Func& f, T1& a1, const T2& a2, PData& pdat);

};


struct Func1Base
    {
    Func1Base() { }
    virtual ~Func1Base() { }

    REGISTER(ITResult virtual operator(),,=0)

    template <typename T>
    ITResult
    operator()(T& t)
        {
        throw ITError("Operation not defined for ITData subtype. [Func1Base]");
        return ITResult();
        }
    };

template <typename Derived>
class Func1Dispatch : public Func1Base
    {
    public:
    Func1Dispatch() { }
    virtual ~Func1Dispatch() { }

    REGISTER(ITResult operator(),,final { return static_cast<Derived*>(this)->applyTo(t); })
    };

template <typename Callable>
class Func1 : public Func1Dispatch<Func1<Callable>>
    {
    Callable& d_;
    PData& pdat_;
    public:
    Func1(Callable& d, PData& pdat) : d_(d), pdat_(pdat) { }
    virtual ~Func1() { }

    template <typename DataType>
    ITResult
    applyTo(DataType& t) { return detail::clone_modify(d_,t,pdat_); }
    };

template <typename Callable>
class Func2Mod : public Func1Dispatch<Func2Mod<Callable>>
    {
    Callable& d_;
    const ITData& arg2_;
    PData& pdat1_;
    public:

    Func2Mod(Callable& d, const ITData& arg2, PData& pdat1) : d_(d), arg2_(arg2), pdat1_(pdat1) { }
    virtual ~Func2Mod() { }

    template<typename DataType>
    ITResult
    applyTo(DataType& arg1);
    };

//////////////////
//////////////////


struct ConstFunc1Base
    {
    ConstFunc1Base() { }
    virtual ~ConstFunc1Base() { }

    REGISTER(ITResult virtual operator(),const,=0)

    template <typename T>
    ITResult
    operator()(const T& t)
        {
        throw ITError("Operation not defined for ITData subtype. [ConstFunc1Base]");
        return ITResult();
        }
    };

template <typename Derived>
class ConstFunc1Dispatch : public ConstFunc1Base
    {
    public:
    ConstFunc1Dispatch() { }
    virtual ~ConstFunc1Dispatch() { }

    REGISTER(ITResult operator(),const,final { return static_cast<Derived*>(this)->applyTo(t); })
    };


template <typename Callable>
class ConstFunc1 : public ConstFunc1Dispatch<ConstFunc1<Callable>>
    {
    Callable& d_;
    public:
    ConstFunc1(Callable& d) : d_(d) { }
    virtual ~ConstFunc1() { }

    template <typename DataType>
    ITResult 
    applyTo(const DataType& t) { return detail::call<ITResult>(d_,t); }
    };

template <typename Callable>
class Func2 : public ConstFunc1Dispatch<Func2<Callable>>
    {
    Callable& d_;
    const ITData& arg2_;
    public:

    Func2(Callable& d, const ITData& arg2) : d_(d), arg2_(arg2) { }
    virtual ~Func2() { }

    template<typename DataType>
    ITResult
    applyTo(const DataType& arg1);
    };


//////////////////
//////////////////


class ITData
    {
    public:

    ITData() { }
    virtual ~ITData() { }

    NewData virtual
    clone() const = 0;

    ITResult virtual
    plugInto(ConstFunc1Base& f) const = 0;

    ITResult virtual
    plugInto(Func1Base& f) = 0;
    };

template <class Derived>
struct ITDispatch : public ITData
    {
    ITDispatch() { }
    virtual ~ITDispatch() { }

    private:
    
    NewData
    clone() const final 
        { 
        auto pdt = static_cast<const Derived*>(this);
        return std::make_unique<Derived>(*pdt);
        }

    ITResult
    plugInto(ConstFunc1Base& f) const final
        {
        const Derived& dt = *(static_cast<const Derived*>(this));
        return f(dt);
        }
        
    ITResult
    plugInto(Func1Base& f) final
        {
        Derived& dt = *(static_cast<Derived*>(this));
        return f(dt);
        }
    };

//
// Implementations
//

namespace detail {

template<typename Func, typename T>
auto
clone_modify_impl(Func& f, T& a, PData& pdat,int) -> decltype(f(static_cast<const T&>(a)))
    {
    println("Not calling solo");
    const T& ca = a;
    return f(ca);
    }

template<typename Func, typename T>
ITResult
clone_modify_impl(Func& f, T& a, PData& pdat,long)
    {
    println("Calling solo");
    T *pa = &a;
    if(!pdat.unique()) 
        {
        pdat = pdat->clone();
        pa = static_cast<T*>(pdat.get());
        }
    return detail::call<ITResult>(f,*pa);
    }

template<typename Func, typename T>
ITResult
clone_modify(Func& f, T& a, PData& pdat)
    {
    return clone_modify_impl(f,a,pdat,0);
    }

/////////////////////

template<typename Func, typename T1, typename T2>
auto
clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,int) -> decltype(f(static_cast<const T1&>(a1)))
    {
    println("Not calling solo (2 arg version)");
    const T1& ca1 = a1;
    return f(ca1,a2);
    }

template<typename Func, typename T1, typename T2>
ITResult
clone_modify_impl(Func& f, T1& a1, const T2& a2, PData& pdat,long)
    {
    println("Calling solo (2 arg version)");
    T1 *pa1 = &a1;
    if(!pdat.unique()) 
        {
        pdat = pdat->clone();
        pa1 = static_cast<T1*>(pdat.get());
        }
    return detail::call<ITResult>(f,*pa1,a2);
    }

template<typename Func, typename T1, typename T2>
ITResult
clone_modify(Func& f, T1& a1, const T2& a2, PData& pdat)
    {
    return clone_modify_impl(f,a1,a2,pdat,0);
    }

};

template<typename Callable>
template<typename DataType>
ITResult Func2<Callable>::
applyTo(const DataType& arg1)
    {
    auto C = [this,&arg1](const auto& a2) { return detail::call<ITResult>(this->d_,arg1,a2); };
    auto f1 = ConstFunc1<decltype(C)>(C);
    return arg2_.plugInto(f1);
    }

template<typename Callable>
template<typename DataType>
ITResult Func2Mod<Callable>::
applyTo(DataType& arg1)
    {
    auto C = [this,&arg1](const auto& a2) { return detail::call<ITResult>(this->d_,arg1,a2); };
    auto f1 = ConstFunc1<decltype(C)>(C);
    return arg2_.plugInto(f1);
    }


//
// applyFunc methods
//

template<typename F>
F
applyFunc(const ITData& arg,
          F f = F())
    {
    ConstFunc1<F> cf1(f);
    arg.plugInto(cf1);
    return f;
    }

template<typename F>
F
applyFunc(const CPData& arg,
          F f = F())
    {
    return applyFunc(*arg,f);
    }

template<typename F>
F
applyFunc(PData& arg,
          F f = F())
    {
    Func1<F> f1(f,arg);
    auto res = arg->plugInto(f1);
    res.update(arg);
    return f;
    }

template<typename F>
F
applyFunc(const PData& arg1,
          const PData& arg2,
          F f = F())
    {
    Func2<F> f2(f,*arg1);
    arg2->plugInto(f2);
    return f;
    }

template<typename F>
F 
applyFunc(PData& arg1,
          const PData& arg2,
          F f = F())
    {
    Func2Mod<F> f2m(f,*arg2,arg1);
    auto res = arg1->plugInto(f2m);
    res.update(arg1,arg2);
    return f;
    }

}; //namespace itensor

#endif

