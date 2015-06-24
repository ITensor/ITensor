//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DOTASK_H
#define __ITENSOR_DOTASK_H

#include <cassert>
#include "itensor/itdata/itdata.h"
//#define REGISTER_ITDATA_HEADER_FILES
//#include "itensor/itdata/storage_types.h"

namespace itensor {

///////////////////

namespace detail {

//Some definitions to help simplify
//template overload selection
//(taken from blog post by R. Martinho Fernandes)

template<unsigned I>
struct choice : choice<I+1> { };

template<>
struct choice<10> { };

struct select_overload : choice<1> { };

struct otherwise{ otherwise(...){} };


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
applyFunc_impl(otherwise, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m)
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
template<typename Return>
struct OneArg
    {
    template<typename Task, typename D>
    Return
    call(Task& t, const D& d, ManageStore& m);

    template<typename Task, typename D>
    Return
    call(Task& t, D& d, ManageStore& m);
    };

template<typename Return>
struct TwoArgs
    {
    template<typename Task, typename D>
    Return
    call(Task& t, const D& d, ManageStore& m);

    template<typename Task, typename D>
    Return
    call(Task& t, D& d, ManageStore& m);
    };

template<typename Derived, typename List>
struct FuncT : FuncT<Derived,typename List::Next>
    {
    using T = typename List::Type;
    using FuncT<Derived,typename List::Next>::applyTo;

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
getReturnHelperImpl(otherwise, RegisterTask<NArgs,Task,Return>& R)
    {
    return std::move(R.task_);
    }
template<typename NArgs, typename Task, typename Return>
auto
getReturnHelperImpl(choice<1>, RegisterTask<NArgs,Task,Return>& R)
    -> typename std::enable_if<std::is_same<typename RegisterTask<NArgs,Task,Return>::return_type,Return>::value,Return>::type
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

    template<typename... MSArgs>
    RegisterTask(task_type&& t,
                 MSArgs&&... msargs)
      : task_(std::move(t)),
        m_(std::forward<MSArgs>(msargs)...)
        { }

    RegisterTask(RegisterTask&& o)
      : task_(std::move(o.task_)), 
        m_(std::move(o.m_)),
        ret_(std::move(o.ret_))
        { }

    virtual ~RegisterTask() { }

    return_type
    getReturn() { return getReturnHelper(*this); }

    template<typename D>
    void
    applyToImpl(const D& d)
        {
        //bool has_evaluate = tryEvaluate(d,m);
        //if(has_evaluate)
        //    {
        //    }
        //else
        ret_ = std::move(NArgs{}.call(task_,d,m_));
        }
    template<typename D>
    void
    applyToImpl(D& d)
        {
        //bool has_evaluate = tryEvaluate(d,m);
        //if(has_evaluate)
        //    {
        //    }
        //else
        ret_ = std::move(NArgs{}.call(task_,d,m_));
        }
    };

template <typename Task, typename D1, typename Return>
struct CallWrap : public FuncT<CallWrap<Task,D1,Return>,StorageTypes>
    {
    CallWrap(Task& t, D1& arg1, ManageStore& m) 
        : t_(t), arg1_(arg1), parg1_(nullptr), m_(m) { }
    CallWrap(Task& t, D1& arg1, PData& parg1, ManageStore& m) 
        : t_(t), arg1_(arg1), parg1_(&parg1), m_(m) { }

    Return
    getReturn() { return std::move(ret_); }


    template<typename D2>
    void
    applyToImpl(const D2& d2);
    
    private:

    Task& t_;
    D1& arg1_;
    PData* parg1_;
    Return ret_;
    ManageStore& m_;

    };

//
// Implementations
//

//If ActualRet!=void this gets called
template<typename Ret, typename ActualRet, typename... VArgs>
struct FixRet
    {
    Ret
    operator()(VArgs&&... vargs)
        {
        return doTask(std::forward<VArgs>(vargs)...);
        }
    };
//If ActualRet==void this gets called
template<typename Ret, typename... VArgs>
struct FixRet<Ret,void,VArgs...>
    {
    NoneType
    operator()(VArgs&&... vargs)
        {
        doTask(std::forward<VArgs>(vargs)...);
        return NoneType{};
        }
    };

template<typename Ret>
using RetOrNone = std::conditional_t<std::is_void<Ret>::value,NoneType,Ret>;

template<typename Ret, typename ActualRet>
using RetIfExists = std::conditional_t<std::is_void<ActualRet>::value,Ret,Ret>;

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

template <typename Ret, typename Task, typename Storage>
Ret
callDoTask_Impl(choice<3>, Task& t, Storage& s, ManageStore& m)
    {
    Error("1 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
template <typename Ret,typename Task, typename Storage>
auto 
callDoTask_Impl(choice<2>, Task& t, Storage& s, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,s))>
    {
    return FixRet<Ret,decltype(doTask(t,s)),Task&,Storage&>()(t,s);
    }
//Template substitution will fail (not a compile-time error)
//if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret, typename Task, typename Storage>
auto 
callDoTask_Impl(choice<1>, Task& t, Storage& s, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,s,m))>
    {
    return FixRet<Ret,decltype(doTask(t,s,m)),Task&,Storage&,ManageStore&>()(t,s,m);
    }
//This version of callDoTask attempts "return doTask(t,s,mp);"
//- If doTask(t,s,mp) not defined, tries calling doTask(t,s)
//  then aborts if that is not defined
//  (this converts what would otherwise be 
//  a compile-time error into a run-time error)
//- If doTask(t,s,mp) or doTask(t,s) defined but 
//  either returns void, returns NoneType{} instead
template<typename Ret, typename Task, typename Storage>
Ret
callDoTask(Task& t, Storage& s, ManageStore& m)
    {
    return callDoTask_Impl<Ret,Task,Storage>(select_overload{},t,s,m);
    }

/////////////////////////////////////////////////////

template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_impl(choice<3>,Task& t, const D& cd, ManageStore& m)
    {
    auto* pd = m.modifyData(cd);
    return callDoTask<Ret>(t,*pd,m);
    }
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_impl(choice<2>,Task& t, const D& cd, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,cd))>
    {
    return FixRet<Ret,decltype(doTask(t,cd)),Task&,const D&>()(t,cd);
    }
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_impl(choice<1>,Task& t, const D& cd, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,cd,m))>
    {
    return FixRet<Ret,decltype(doTask(t,cd,m)),Task&,const D&,ManageStore&>()(t,cd,m);
    }
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask(Task& t, const D& cd, ManageStore& m)
    {
    return cloneDoTask_impl<Ret>(select_overload{},t,cd,m);
    }

/////////////////////

//(4) Least preferred version which always compiles
template <typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask_Impl(choice<3>, Task& t, D1& d1, const D2& d2, ManageStore& m)
    {
    Error("2 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(choice<2>, Task& t, D1& d1, const D2& d2, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,d1,d2))>
    {
    return FixRet<Ret,decltype(doTask(t,d1,d2)),Task&,D1&,const D2&>()(t,d1,d2);
    }
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(choice<1>, Task& t, D1& d1, const D2& d2, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,d1,d2,m))>
    {
    return FixRet<Ret,decltype(doTask(t,d1,d2,m)),Task&,D1&,const D2&,ManageStore&>()(t,d1,d2,m);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask(Task& t, D1& d1, const D2& d2, ManageStore& m)
    {
    return callDoTask_Impl<Ret,Task,D1,D2>(select_overload{},t,d1,d2,m);
    }

/////////////////////

template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_impl(choice<3>, Task& t, const D1& cd1, const D2& d2, ManageStore& m)
    {
    auto* pd1 = m.modifyData(cd1);
    return callDoTask<Ret>(t,*pd1,d2,m);
    }
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_impl(choice<2>, Task& t, const D1& cd1, const D2& d2, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,cd1,d2))>
    {
    return FixRet<Ret,decltype(doTask(t,cd1,d2)),Task&,const D1&,const D2&>()(t,cd1,d2);
    }
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_impl(choice<1>, Task& t, const D1& cd1, const D2& d2, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,cd1,d2,m))>
    {
    return FixRet<Ret,decltype(doTask(t,cd1,d2,m)),Task&,const D1&,const D2&,ManageStore&>()(t,cd1,d2,m);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask(Task& t, const D1& d1, const D2& d2, ManageStore& m)
    {
    return cloneDoTask_impl<Ret,Task,D1,D2>(select_overload{},t,d1,d2,m);
    }

void inline
check(const CPData& p)
    {
    if(!p) Error("doTask called on unallocated store pointer");
    }


template<typename Ret>
template<typename Task, typename D>
Ret OneArg<Ret>::
call(Task& t, const D& d, ManageStore& m)
    {
    return detail::callDoTask<Ret>(t,d,m);
    }

template<typename Ret>
template<typename Task, typename D>
Ret OneArg<Ret>::
call(Task& t, D& d, ManageStore& m)
    {
    return detail::cloneDoTask<Ret>(t,d,m);
    }

template<typename Ret>
template<typename Task, typename D>
Ret TwoArgs<Ret>::
call(Task& t, const D& d, ManageStore& m)
    {
    assert(m.hasArg2());
    CallWrap<Task,const D,Ret> w(t,d,m);
    m.arg2().plugInto(w);
    return w.getReturn();
    }

template<typename Ret>
template<typename Task, typename D>
Ret TwoArgs<Ret>::
call(Task& t, D& d, ManageStore& m)
    {
    assert(m.hasPArg1());
    assert(m.hasArg2());
    CallWrap<Task,D,Ret> w(t,d,m.parg1(),m);
    m.arg2().plugInto(w);
    return w.getReturn();
    }

template <typename Task, typename D1, typename Return>
template<typename D2>
void CallWrap<Task,D1,Return>::
applyToImpl(const D2& d2)
    { 
    if(parg1_) ret_ = detail::cloneDoTask<Return>(t_,arg1_,d2,m_); 
    else       ret_ = detail::callDoTask<Return>(t_,arg1_,d2,m_); 
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
       const CPData& arg)
    {
#ifdef DEBUG
    detail::check(arg);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    detail::RegisterTask<detail::OneArg<Ret>,Task,Ret> r(std::forward<Task>(t));
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
    detail::RegisterTask<detail::OneArg<Ret>,Task,Ret> r(std::forward<Task>(t),&arg);
    arg->plugInto(r);
    return std::move(r.getReturn());
    }

template<typename Task>
auto
doTask(Task&& t,
       const CPData& arg1,
       const CPData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    detail::RegisterTask<detail::TwoArgs<Ret>,Task,Ret> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getReturn());
    }

template<typename Task>
auto
doTask(Task&& t,
       PData& arg1,
       const CPData& arg2)
    {
#ifdef DEBUG
    detail::check(arg1);
    detail::check(arg2);
#endif
    using Ret = ReturnType<Task,StorageTypes>;
    detail::RegisterTask<detail::TwoArgs<Ret>,Task,Ret> r(std::forward<Task>(t),&arg1,&arg2);
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

