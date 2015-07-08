//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DOTASK_H
#define __ITENSOR_DOTASK_H

#include <cassert>
#include "itensor/itdata/itdata.h"
#include "itensor/util/call_if.h"
#include "itensor/util/print.h"

namespace itensor {

///////////////////

namespace detail {

//Some definitions to help simplify
//template overload selection
//(credit to R. Martinho Fernandes)

template<unsigned I>
struct choice : choice<I+1> { constexpr choice(){} };

template<>
struct choice<10> { constexpr choice(){} };

struct select_overload : choice<1> { constexpr select_overload(){} };


template<typename F,typename Ret>
struct ApplyFunc
    { 
    using function_type = F;
    F& f;
    Ret r;
    ApplyFunc(F&& f_) : f(f_) { }
    template<typename S>
    void
    operator()(S& s) 
        { 
        r = f(s); 
        }
    operator Ret() const { return std::move(r); }
    };

template<typename F>
struct ApplyFunc<F,void>
    { 
    using function_type = F;
    F& f;
    ApplyFunc(F&& f_) : f(f_) { }
    template<typename S>
    void
    operator()(S& s) { f(s); }
    };

template<typename F, typename R, typename Storage>
void
applyFunc_impl(choice<3>, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m)
    {
    throw ITError("applyFunc: function object has no operator() method for storage type");
    }

template<typename F, typename R, typename Storage>
auto
applyFunc_impl(choice<2>, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m)
    -> decltype(A.f(s), void())
    {
    Storage& ncs = m.modifyData();
    A(ncs);
    }

template<typename F, typename R, typename Storage>
auto
applyFunc_impl(choice<1>, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m) 
    -> decltype(A.f(s), void())
    {
    A(s);
    }

} //namespace detail

template<typename F, typename R, typename Storage>
void
doTask(detail::ApplyFunc<F,R>& A, const Storage& s, ManageStore& m) 
    { 
    detail::applyFunc_impl(detail::select_overload{},A,s,m);
    }


namespace detail {

//OneArg and TwoArgs are "policy classes" 
//for customizing implementation of RegisterTask
template<class PType>
struct OneArg
    {
    using ptype = PType;

    //template<typename RT, typename Task, typename D, typename Return>
    //void
    //call(RT& rt, Task& t, const D& d, ManageStore& m, Return& ret);

    template<typename RT, typename Task, typename D, typename Return>
    void
    call(RT&, Task& t, D& d, ManageStore& m, Return& ret);
    };

template<class PType1, class PType2>
struct TwoArgs
    {
    using ptype1 = PType1;
    using ptype2 = PType2;

    //template<typename RT, typename Task, typename D, typename Return>
    //void
    //call(RT& rt, Task& t, const D& d, ManageStore& m, Return& ret);

    template<typename RT, typename Task, typename D, typename Return>
    void
    call(RT& rt, Task& t, D& d, ManageStore& m, Return& ret);
    };

template<typename Derived, typename TList>
struct FuncT : FuncT<Derived,popFront<TList>>
    {
    using T = frontType<TList>;
    using FuncT<Derived,popFront<TList>>::applyTo;

    void
    applyTo(const T& t) final
        {
        auto* pd = static_cast<Derived*>(this);
        pd->applyToImpl(t);
        }

    void
    applyTo(T& t) final
        {
        auto* pd = static_cast<Derived*>(this);
        pd->applyToImpl(t);
        }
    };
template<typename Derived>
struct FuncT<Derived,TypeList<>> : FuncBase
    { };

template <typename NArgs, typename Task, typename Return>
class RegisterTask;

template<typename NArgs, typename Task, typename Return>
Task
getReturnHelperImpl(choice<2>, RegisterTask<NArgs,Task,Return>& R)
    {
    return std::move(R.task_);
    }
template<typename NArgs, typename Task, typename Return>
auto
getReturnHelperImpl(choice<1>, RegisterTask<NArgs,Task,Return>& R)
    -> std::enable_if_t<std::is_same<typename RegisterTask<NArgs,Task,Return>::return_type,Return>::value,Return>
    {
    return std::move(R.ret_);
    }
template<typename NArgs, typename Task, typename Return>
auto
getReturnHelper(RegisterTask<NArgs,Task,Return>& R)
    {
    return getReturnHelperImpl(select_overload{},R);
    }


template <typename NArgs, typename Task, typename Return>
class RegisterTask : public FuncT<RegisterTask<NArgs,Task,Return>,StorageTypes>
    {
    public:
    using task_type = std::remove_reference_t<Task>;
    using return_type = std::conditional_t<std::is_same<Return,NoneType>::value,
                                           task_type,
                                           Return>;
    task_type task_;
    ManageStore m_;
    Return ret_;

    RegisterTask(task_type&& t)
      : task_(std::move(t))
        { }

    RegisterTask(task_type&& t,
                 ManageStore&& m)
      : task_(std::move(t)),
        m_(std::move(m))
        { }

    virtual ~RegisterTask() { }

    return_type
    getReturn() { return getReturnHelper(*this); }

    template<typename D>
    void
    applyToImpl(D& d);
    };

template <class RT, typename Task, typename D1, typename Return, class PType1, class PType2>
class CallWrap : public FuncT<CallWrap<RT,Task,D1,Return,PType1,PType2>,StorageTypes>
    {
    RT& rt_;
    Task& t_;
    D1& d1_;
    ManageStore& m_;
    Return& ret_;
    public:

    using ptype1 = PType1;
    using ptype2 = PType2;

    CallWrap(RT& rt, Task& t, D1& arg1, ManageStore& m, Return& ret) 
        : rt_(rt), t_(t), d1_(arg1), m_(m), ret_(ret) { }

    template<typename D2>
    void
    applyToImpl(D2& d2);
    
    };

//
// Implementations
//

template<typename Ret>
using RetOrNone = std::conditional_t<std::is_void<Ret>::value,NoneType,Ret>;

template<typename T, typename Ret>
using ReturnIfExists = std::conditional_t<std::is_void<T>::value,Ret,Ret>;

template<typename Task, typename Storage>
NoneType
testRetImpl(choice<3>, Task& t, Storage& s, ManageStore& m)
    {
    return NoneType{};
    }
template<typename Task, typename Storage>
auto 
testRetImpl(choice<2>, Task& t, Storage& s, ManageStore& m)
    -> RetOrNone<decltype(doTask(t,s))>
    {
    using Ret = RetOrNone<decltype(doTask(t,s))>;
    return Ret{};
    }
template<typename Task, typename Storage>
auto 
testRetImpl(choice<1>, Task& t, Storage& s, ManageStore& m)
    -> RetOrNone<decltype(doTask(t,s,m))>
    {
    using Ret = RetOrNone<decltype(doTask(t,s,m))>;
    return Ret{};
    }
template<typename Task, typename Storage>
struct TestRet
    {
    Task& t;
    Storage& s;
    ManageStore& m;
    TestRet(Task& t_, Storage& s_, ManageStore& m_)
        : t(t_), s(s_), m(m_) 
        { }
    auto
    operator()()
        {
        return testRetImpl(select_overload{},t,s,m);
        }
    };

/////////////

//If ActualRet!=void this gets called
template<typename Ret, typename ActualRet>
class FixRet
    {
    Ret& ret_;
    public:

    FixRet(Ret& ret) : ret_(ret) { }

    template<typename... VArgs>
    void
    operator()(VArgs&&... vargs)
        {
        ret_ = doTask(std::forward<VArgs>(vargs)...);
        }
    };
//If ActualRet==void this gets called
template<typename Ret>
class FixRet<Ret,void>
    {
    public:

    FixRet(Ret& ret) { }

    template<typename... VArgs>
    void
    operator()(VArgs&&... vargs)
        {
        doTask(std::forward<VArgs>(vargs)...);
        }
    };

/////////////

template <typename Ret, typename Task, typename Storage>
void
callDoTask_Impl(choice<3>, Task& t, Storage& s, ManageStore& m, Ret& ret)
    {
    throw ITError("1 parameter doTask not defined for specified task or data type [1]");
    }
template <typename Ret,typename Task, typename Storage>
auto 
callDoTask_Impl(choice<2>, Task& t, Storage& s, ManageStore& m, Ret& ret)
    -> ReturnIfExists<decltype(doTask(t,s)),void>
    {
    FixRet<Ret,decltype(doTask(t,s))>{ret}(t,s);
    }
template <typename Ret, typename Task, typename Storage>
auto 
callDoTask_Impl(choice<1>, Task& t, Storage& s, ManageStore& m, Ret& ret)
    -> ReturnIfExists<decltype(doTask(t,s,m)),void>
    {
    FixRet<Ret,decltype(doTask(t,s,m))>{ret}(t,s,m);
    }
template<typename Return, typename Storage, typename Task>
void
callDoTask(Task& t, Storage& s, ManageStore& m, Return& ret)
    {
    callDoTask_Impl<Return,Task,Storage>(select_overload{},t,s,m,ret);
    }

/////////////////////////////////////////////////////

template <typename Ret, typename Task, typename D1, typename D2>
void
callDoTask_Impl(choice<3>, Task& t, D1& d1, const D2& d2, ManageStore& m, Ret& ret)
    {
    throw ITError("2 parameter doTask not defined for specified task or data type [2]");
    }
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(choice<2>, Task& t, D1& d1, const D2& d2, ManageStore& m, Ret& ret)
    -> ReturnIfExists<decltype(doTask(t,d1,d2)),void>
    {
    FixRet<Ret,decltype(doTask(t,d1,d2))>{ret}(t,d1,d2);
    }
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(choice<1>, Task& t, D1& d1, const D2& d2, ManageStore& m, Ret& ret)
    -> ReturnIfExists<decltype(doTask(t,d1,d2,m)),void>
    {
    FixRet<Ret,decltype(doTask(t,d1,d2,m))>{ret}(t,d1,d2,m);
    }
template<typename Ret, typename Task, typename D1, typename D2>
void
callDoTask(Task& t, D1& d1, const D2& d2, ManageStore& m, Ret& ret)
    {
    callDoTask_Impl<Ret,Task,D1,D2>(select_overload{},t,d1,d2,m,ret);
    }

/////////////////////

template<typename Task, typename D>
std::false_type
testDTImpl(choice<3>, Task& t, D& d, ManageStore& m)
    {
    return std::false_type{};
    }
template<typename Task, typename D>
auto 
testDTImpl(choice<2>, Task& t, D& d, ManageStore& m)
    -> ReturnIfExists<decltype(doTask(t,d)),std::true_type>
    {
    return std::true_type{};
    }
template<typename Task, typename D>
auto 
testDTImpl(choice<1>, Task& t, D& d, ManageStore& m)
    -> ReturnIfExists<decltype(doTask(t,d,m)),std::true_type>
    {
    return std::true_type{};
    }
template<typename Task, typename Storage>
struct HasDTHelper
    {
    Task* t;
    Storage* s;
    ManageStore* m;
    auto operator()() { return testDTImpl(select_overload{},*t,*s,*m); }
    };
template<typename Task, typename Storage>
struct HasConstDoTask
    {
    using ResultType = std::result_of_t<HasDTHelper<Task,const Storage>()>;
    bool constexpr static
    result() { return ResultType{}; }
    };
template<typename Task, typename Storage>
struct HasNonConstDoTask
    {
    using ResultType = std::result_of_t<HasDTHelper<Task,std::remove_const_t<Storage>>()>;
    using CResultType = std::result_of_t<HasDTHelper<Task,const Storage>()>;
    bool constexpr static
    result() { return ResultType{} && (not CResultType{}); }
    };
template<typename Task, typename Storage>
struct HasDoTask
    {
    using CResultType = std::result_of_t<HasDTHelper<Task,const Storage>()>;
    using NCResultType = std::result_of_t<HasDTHelper<Task,std::remove_const_t<Storage>>()>;
    using ResultType = std::conditional_t<CResultType::value,
                                          CResultType,
                                          NCResultType>;
    bool constexpr static
    result() { return ResultType{}; }
    };

/////////////////////

template<typename Task, typename D1, typename D2>
std::false_type
testDTImpl(choice<3>, Task& t, D1& d1, const D2& d2, ManageStore& m)
    {
    return std::false_type{};
    }
template<typename Task, typename D1, typename D2>
auto 
testDTImpl(choice<2>, Task& t, D1& d1, const D2& d2, ManageStore& m)
    -> ReturnIfExists<decltype(doTask(t,d1,d2)),std::true_type>
    {
    return std::true_type{};
    }
template<typename Task, typename D1, typename D2>
auto 
testDTImpl(choice<1>, Task& t, D1& d1, const D2& d2, ManageStore& m)
    -> ReturnIfExists<decltype(doTask(t,d1,d2,m)),std::true_type>
    {
    return std::true_type{};
    }
template<typename Task, typename D1, typename D2>
struct HasDTHelper2Arg
    {
    Task* t;
    D1* d1;
    D2* d2;
    ManageStore* m;
    auto operator()() { return testDTImpl(select_overload{},*t,*d1,*d2,*m); }
    };
template<typename Task, typename D1, typename D2>
struct HasConstDoTask2Arg
    {
    using ResultType = std::result_of_t<HasDTHelper2Arg<Task,const D1, const D2>()>;
    bool constexpr static
    result() { return ResultType{}; }
    };
template<typename Task, typename D1, typename D2>
struct HasNonConstDoTask2Arg
    {
    using ResultType = std::result_of_t<HasDTHelper2Arg<Task,std::remove_const_t<D1>,D2>()>;
    using CResultType = std::result_of_t<HasDTHelper2Arg<Task,const D1,D2>()>;
    bool constexpr static
    result() { return ResultType{} && (not CResultType{}); }
    };
template<typename Task, typename D1, typename D2>
struct HasDoTask2Arg
    {
    using CResultType = std::result_of_t<HasDTHelper2Arg<Task,const D1,D2>()>;
    using NCResultType = std::result_of_t<HasDTHelper2Arg<Task,std::remove_const_t<D1>,D2>()>;
    using ResultType = std::conditional_t<CResultType::value,
                                          CResultType,
                                          NCResultType>;
    bool constexpr static
    result() { return ResultType{}; }
    };

/////////////////////

template<typename D>
auto 
testEvalImpl(choice<2>, D& d)
    {
    return std::false_type{};
    }
template<typename D>
auto 
testEvalImpl(choice<1>, D& d)
    -> ReturnIfExists<decltype(evaluate(d)),std::true_type>
    {
    return std::true_type{};
    }
template<typename Storage>
struct HasEvaluate
    {
    struct Test 
        {
        Storage* s;
        auto operator()() { return testEvalImpl(select_overload{},*s); }
        };
    bool constexpr static
    result() { return std::result_of_t<Test()>{}; }
    };

/////////////

template<typename D>
PData
callEvaluateImpl(choice<2>, D& d)
    {
    throw std::runtime_error("No doTask overload found for task/storage type");
    return PData{};
    }
template<typename D>
auto
callEvaluateImpl(choice<1>, D& d)
    -> ReturnIfExists<decltype(evaluate(d)),PData>
    {
    return evaluate(d);
    }
template<typename D>
PData
callEvaluate(D& d)
    {
    return callEvaluateImpl(select_overload{},d);
    }

/////////////

template<typename D>
bool constexpr
checkHasResultImpl(choice<2>, const D& d)
    {
    return true;
    }
template<typename D>
auto
checkHasResultImpl(choice<1>, const D& d)
    -> ReturnIfExists<decltype(hasResult(d)),bool>
    {
    return hasResult(d);
    }
template<typename D>
bool constexpr
checkHasResult(const D& d)
    {
    return checkHasResultImpl(select_overload{},d);
    }

/////////////


void inline
check(const PData& p)
    {
    if(!p) Error("doTask called on unallocated store pointer");
    }
void inline
check(const CPData& p)
    {
    if(!p) Error("doTask called on unallocated store pointer");
    }

/////////////

template <typename NArgs, typename Task, typename Return>
template<typename D>
void RegisterTask<NArgs,Task,Return>::
applyToImpl(D& d)
    {
    constexpr bool isLazy = HasEvaluate<D>::result();
    if(isLazy)
        {
        if(checkHasResult(d))
            {
            println("Lazy storage has result, swapping");
            m_.parg1() = callEvaluate(d);
            m_.parg1()->plugInto(*this);
            return;
            }
        }
    NArgs{}.call(*this,task_,d,m_,ret_);
    }

template<class PType>
template<typename RT, typename Task, typename D, typename Return>
void OneArg<PType>::
call(RT& rt, Task& t, D& d, ManageStore& m, Return& ret)
    {
    constexpr bool NCData = std::is_same<PType,PData>::value;
    constexpr bool hasConstDT = HasConstDoTask<Task,D>::result();
    constexpr bool hasNCDT = HasNonConstDoTask<Task,D>::result();
    constexpr bool isLazy = HasEvaluate<D>::result();
    if(hasConstDT)
        {
        const auto& cd = d;
        detail::callDoTask(t,cd,m,ret);
        }
    else if(NCData && hasNCDT)
        {
        auto* pd = m.modifyData(d);
        detail::callDoTask(t,*pd,m,ret);
        }
    else if(isLazy)
        {
        m.parg1() = callEvaluate(d);
        m.parg1()->plugInto(rt);
        }
    else
        {
        throw ITError("doTask not defined for task/storage type [4]");
        }
    }

template<class PType1, class PType2>
template<typename RT, typename Task, typename D, typename Return>
void TwoArgs<PType1,PType2>::
call(RT& rt, Task& t, D& d, ManageStore& m, Return& ret)
    {
    CallWrap<RT,Task,D,Return,PType1,PType2> w(rt,t,d,m,ret);
    m.parg2()->plugInto(w);
    }

template <class RT, typename Task, typename D1, typename Return, class PType1, class PType2>
template<typename D2>
void CallWrap<RT,Task,D1,Return,PType1,PType2>::
applyToImpl(D2& d2)
    { 
    constexpr bool isLazy2 = HasEvaluate<D2>::result();
    if(isLazy2)
        {
        if(checkHasResult(d2))
            {
            println("Lazy storage has result, swapping");
            m_.parg2() = callEvaluate(d2);
            m_.parg2()->plugInto(*this);
            return;
            }
        }

    constexpr bool NCData1 = std::is_same<PType1,PData>::value;
    constexpr bool hasCDT = HasConstDoTask2Arg<Task,D1,D2>::result();
    constexpr bool hasNCDT = HasNonConstDoTask2Arg<Task,D1,D2>::result();
    constexpr bool isLazy1 = HasEvaluate<D1>::result();
    if(hasCDT)
        {
        const auto& cd1 = d1_;
        callDoTask(t_,cd1,d2,m_,ret_);
        }
    else if(NCData1 && hasNCDT)
        {
        auto* pd1 = m_.modifyData(d1_);
        callDoTask(t_,*pd1,d2,m_,ret_);
        }
    else if(isLazy1 && isLazy2)
        {
        m_.parg1() = callEvaluate(d1_);
        m_.parg2() = callEvaluate(d2);
        m_.parg1()->plugInto(rt_);
        }
    else if(isLazy1)
        {
        constexpr bool hasCPDT = HasConstDoTask2Arg<Task,D1,const PData>::result();
        constexpr bool hasNCPDT = HasNonConstDoTask2Arg<Task,D1,const PData>::result();
        if(hasCPDT)
            {
            const auto& cd1 = d1_;
            callDoTask(t_,cd1,m_.parg2(),m_,ret_);
            }
        else if(NCData1 && hasNCPDT)
            {
            auto* pd1 = m_.modifyData(d1_);
            callDoTask(t_,*pd1,m_.parg2(),m_,ret_);
            }
        else
            {
            m_.parg1() = callEvaluate(d1_);
            m_.parg1()->plugInto(rt_);
            }
        }
    else if(isLazy2)
        {
        constexpr bool hasPDT = HasConstDoTask2Arg<Task,const PData,D2>::result();
        if(hasPDT)
            {
            callDoTask(t_,m_.parg1(),d2,m_,ret_);
            }
        else
            {
            m_.parg2() = callEvaluate(d2);
            m_.parg2()->plugInto(*this);
            }
        }
    else
        {
        throw ITError("doTask not defined for task/storage type [6]");
        }
    }

template<typename Task, typename TList>
struct GetRType : GetRType<Task,popFront<TList>>
    {
    using Test = std::result_of_t<TestRet<Task,frontType<TList>>()>;
    using Parent = GetRType<Task,popFront<TList>>;
    using RType = std::conditional_t<not std::is_same<Test,NoneType>::value,
                                     Test,
                                     typename Parent::RType>;
    };

template<typename Task>
struct GetRType<Task,TypeList<>>
    {
    using RType = NoneType;
    };

} //namespace detail

template<typename Task, typename TList>
using ReturnType = typename detail::GetRType<std::remove_reference_t<Task>,TList>::RType;

//////
////// doTask methods
//////

template<typename T, typename... VArgs>
PData
newITData(VArgs&&... vargs)
    {
    static_assert(containsType<StorageTypes,T>{},"Data type not in list of registered storage types");
    return std::make_shared<ITWrap<T>>(std::forward<VArgs>(vargs)...);
    }

template<typename Task>
auto
doTask(Task&& t,
       CPData arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    ManageStore m(&(arg.p));
    detail::RegisterTask<detail::OneArg<CPData>,Task,Ret> r(std::forward<Task>(t),std::move(m));
    arg->plugInto(r);
    return std::move(r.getReturn());
    }

template<typename Task>
auto
doTask(Task&& t,
       PData& arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    ManageStore m(&arg);
    detail::RegisterTask<detail::OneArg<PData>,Task,Ret> r(std::forward<Task>(t),std::move(m));
    arg->plugInto(r);
    return std::move(r.getReturn());
    }

template<typename Task>
auto
doTask(Task&& t,
       CPData arg1,
       CPData arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    ManageStore m(&(arg1.p),&(arg2.p));
    detail::RegisterTask<detail::TwoArgs<CPData,CPData>,Task,Ret> r(std::forward<Task>(t),std::move(m));
    arg1->plugInto(r);
    return std::move(r.getReturn());
    }

template<typename Task>
auto
doTask(Task&& t,
       PData& arg1,
       CPData arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    ManageStore m(&arg1,&(arg2.p));
    detail::RegisterTask<detail::TwoArgs<PData,CPData>,Task,Ret> r(std::forward<Task>(t),std::move(m));
    arg1->plugInto(r);
    return std::move(r.getReturn());
    }


template<typename F>
F
applyFunc(F&& f, PData& store)
    {
    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
    return f;
    }

template<typename F>
F
applyFunc(F&& f, const CPData& store)
    {
    doTask(detail::ApplyFunc<F,void>{std::forward<F>(f)},store);
    return f;
    }

template<typename Ret, typename F>
Ret
applyFunc(F&& f, PData& store)
    {
    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
    }

template<typename Ret, typename F>
Ret
applyFunc(F&& f, const CPData& store)
    {
    return doTask(detail::ApplyFunc<F,Ret>{std::forward<F>(f)},store);
    }

} //namespace itensor

#undef REGISTER_ITDATA_HEADER_FILES

#endif

