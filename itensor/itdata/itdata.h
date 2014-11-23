//
// Distributed under the ITensor Library License, Version 1.1
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
    X(Y  ITDiag<Real>       &t) Z;\
    X(Y  ITDiag<Complex>    &t) Z;
///////////////////////////////////

namespace itensor {

///////////////////////////////////
// (2) Forward declarations for all ITData subtypes:
template <typename>
class ITDense;
template <typename>
class ITDiag;
///////////////////////////////////



class ITData;
using PData = std::shared_ptr<ITData>;
using CPData = std::shared_ptr<const ITData>;
using NewData = std::unique_ptr<ITData>;

template<typename DataType, typename... Args>
std::unique_ptr<DataType>
make_newdata(Args&&... args)
    {
    return std::unique_ptr<DataType>(new DataType(std::forward<Args>(args)...));
    }

struct Func1Base
    {
    Func1Base() { }
    virtual ~Func1Base() { }

    REGISTER(NewData virtual operator(),,=0)

    template <typename T>
    NewData
    operator()(T& t)
        {
        throw ITError("Operation not defined for ITData subtype. [Func1Base]");
        return NewData();
        }
    };

template <typename Derived>
class Func1Dispatch : public Func1Base
    {
    public:
    Func1Dispatch() { }
    virtual ~Func1Dispatch() { }

    REGISTER(NewData operator(),,final { return static_cast<Derived*>(this)->applyTo(t); })
    };


struct ConstFunc1Base
    {
    ConstFunc1Base() { }
    virtual ~ConstFunc1Base() { }

    REGISTER(NewData virtual operator(),const,=0)

    template <typename T>
    NewData
    operator()(const T& t)
        {
        throw ITError("Operation not defined for ITData subtype. [ConstFunc1Base]");
        return NewData();
        }
    };

template <typename Derived>
class ConstFunc1Dispatch : public ConstFunc1Base
    {
    public:
    ConstFunc1Dispatch() { }
    virtual ~ConstFunc1Dispatch() { }

    REGISTER(NewData operator(),const,final { return static_cast<Derived*>(this)->applyTo(t); })
    };

template <typename Callable>
class Func1 : public Func1Dispatch<Func1<Callable>>
    {
    Callable& d_;
    public:
    Func1(Callable& d) : d_(d) { }
    virtual ~Func1() { }

    template <typename DataType>
    NewData
    applyTo(DataType& t) { return detail::call<NewData>(d_,t); }
    };

template <typename Callable>
class Func2Mod : public Func1Dispatch<Func2Mod<Callable>>
    {
    Callable& d_;
    const ITData& arg2_;
    public:

    Func2Mod(Callable& d, const ITData& arg2) : d_(d), arg2_(arg2) { }
    virtual ~Func2Mod() { }

    template<typename DataType>
    NewData
    applyTo(DataType& arg1);
    };

template <typename Callable>
class ConstFunc1 : public ConstFunc1Dispatch<ConstFunc1<Callable>>
    {
    Callable& d_;
    public:
    ConstFunc1(Callable& d) : d_(d) { }
    virtual ~ConstFunc1() { }

    template <typename DataType>
    NewData 
    applyTo(const DataType& t) { return detail::call<NewData>(d_,t); }
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
    NewData
    applyTo(const DataType& arg1);
    };

class ITData
    {
    public:

    ITData() { }
    virtual ~ITData() { }

    NewData virtual
    clone() const = 0;

    NewData virtual
    plugInto(ConstFunc1Base& f) const = 0;

    NewData virtual
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

    NewData
    plugInto(ConstFunc1Base& f) const final
        {
        const Derived& dt = *(static_cast<const Derived*>(this));
        return f(dt);
        }
        
    NewData
    plugInto(Func1Base& f) final
        {
        Derived& dt = *(static_cast<Derived*>(this));
        return f(dt);
        }
    };

template<typename Callable>
template<typename DataType>
NewData Func2<Callable>::
applyTo(const DataType& arg1)
    {
    auto C = [this,&arg1](const auto& a2) { return detail::call<NewData>(this->d_,arg1,a2); };
    auto f1 = ConstFunc1<decltype(C)>(C);
    return arg2_.plugInto(f1);
    }

template<typename Callable>
template<typename DataType>
NewData Func2Mod<Callable>::
applyTo(DataType& arg1)
    {
    auto C = [this,&arg1](const auto& a2) { return detail::call<NewData>(this->d_,arg1,a2); };
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

namespace detail {
template<typename F>
NewData
applyFuncImpl(ITData& arg,
              F& f)
    {
    Func1<F> f1(f);
    return arg.plugInto(f1);
    }

};

template<typename F>
F
applyFunc(PData& arg,
          F f = F())
    {
    auto res = detail::applyFuncImpl(*arg,f);
    if(res) arg = std::move(res);
    return f;
    }

template<typename F>
F
applyFunc(NewData& arg,
          F f = F())
    {
    auto res = detail::applyFuncImpl(*arg,f);
    if(res) arg = std::move(res);
    return f;
    }

template<typename F>
F
applyFunc(const ITData& arg1,
          const ITData& arg2,
          F f = F())
    {
    Func2<F> f2(f,arg1);
    arg2.plugInto(f2);
    return f;
    }

template<typename F>
F
applyFunc(const PData& arg1,
          const PData& arg2,
          F f = F())
    {
    return applyFunc(*arg1,*arg2,f);
    }

template<typename F>
F 
applyFunc(PData& arg1,
          const ITData& arg2,
          F f = F())
    {
    Func2Mod<F> f2m(f,arg2);
    auto res = arg1->plugInto(f2m);
    if(res) arg1 = std::move(res);
    return f;
    }

template<typename F>
F 
applyFunc(PData& arg1,
          const PData& arg2,
          F f = F())
    {
    return applyFunc(arg1,*arg2,f);
    }

}; //namespace itensor

#endif

