//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_DOTASK_H
#define __ITENSOR_DOTASK_H

#include "itensor/itdata/itdata.h"
#define REGISTER_ITDATA_HEADER_FILES
#include "itensor/itdata/storage_types.h"

namespace itensor {

///////////////////

namespace detail {

//Some definitions to help document
//the functions below
constexpr const int TryFirst = 0;
using First  = int;
using Second = long;


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
applyFunc_ncimpl(Second, ApplyFunc<F,R>& A, Storage& s)
    {
    throw ITError("applyFunc: function object has no operator() method for storage type");
    }

template<typename F, typename R, typename Storage>
auto
applyFunc_ncimpl(First, ApplyFunc<F,R>& A, Storage& s)
    -> decltype(A.f(s), void())
    {
    A(s);
    }

template<typename F, typename R, typename Storage>
void
applyFunc_constimpl(Second, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m) 
    {
    Storage& ncs = m.modifyData();
    applyFunc_ncimpl(0,A,ncs);
    }

template<typename F, typename R, typename Storage>
auto
applyFunc_constimpl(First, ApplyFunc<F,R>& A, const Storage& s, ManageStore& m) 
    -> decltype(A.f(s), void())
    {
    A(s);
    }

} //namespace detail

template<typename F, typename R, typename Storage>
void
doTask(detail::ApplyFunc<F,R>& A, const Storage& s, ManageStore& m) 
    { 
    detail::applyFunc_constimpl(detail::TryFirst,A,s,m);
    }


///////////////////


struct NoneType { };

//OneArg and TwoArgs are "policy classes" 
//for customizing implementation of RegisterTask
struct OneArg
    {
    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, const D& d, ManageStore& m);

    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, D& d, ManageStore& m);
    };

struct TwoArgs
    {
    template<typename Return, typename Task, typename D>
    Return
    call(Task& t, const D& d, ManageStore& m);

    template<typename Return, typename Task, typename D>
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
    {
    };

template <typename NArgs, typename Task, typename Return>
class RegisterTask : public FuncT<RegisterTask<NArgs,Task,Return>,StorageTypes>
    {
    public:
    using task_type = std::remove_reference_t<Task>;
    using return_type = std::remove_reference_t<Return>;
    private:
    task_type task_;
    ManageStore m_;
    return_type ret_;
    public:

    RegisterTask(task_type&& t) :
        task_(std::move(t))
        { }

    template<typename... MSArgs>
    RegisterTask(task_type&& t,
                 MSArgs&&... msargs) :
        task_(std::move(t)),
        m_(std::forward<MSArgs>(msargs)...)
        { }

    RegisterTask(RegisterTask&& o) :
        task_(std::move(o.task_)), 
        m_(std::move(o.m_)),
        ret_(std::move(o.ret_))
        { }

    virtual ~RegisterTask() { }

    return_type&
    getReturn() { return ret_; }

    task_type&
    getTask() { return task_; }

    template<typename D>
    void
    applyToImpl(const D& d)
        {
        ret_ = std::move(NArgs().template call<return_type>(task_,d,m_));
        }
    template<typename D>
    void
    applyToImpl(D& d)
        {
        ret_ = std::move(NArgs().template call<return_type>(task_,d,m_));
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

namespace detail {

//template<typename Ret, typename... VArgs>
//struct FixRet
//    {
//    Ret
//    operator()(VArgs&&... vargs)
//        {
//        return doTask(std::forward<VArgs>(vargs)...);
//        }
//    };
//template<typename... VArgs>
//struct FixRet<NoneType,VArgs...>
//    {
//    NoneType
//    operator()(VArgs&&... vargs)
//        {
//        doTask(std::forward<VArgs>(vargs)...);
//        return NoneType{};
//        }
//    };

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
struct TestRet
    {
    auto
    operator()()
        {
        return testRetImplMS(TryFirst);
        }
    private:
    Task& t;
    Storage& s;
    ManageStore& m;
    auto 
    testRetImplMS(First)
        -> RetOrNone<decltype(doTask(t,s,m))>
        {
        using ActualRet = RetOrNone<decltype(doTask(t,s,m))>;
        //return FixRet<ActualRet,Task&,Storage&,ManageStore&>()(t,s,m);
        return std::declval<ActualRet>();
        }
    auto
    testRetImplMS(Second)
        {
        return testRetImpl(TryFirst,t,s);
        }
    auto 
    testRetImpl(First)
        -> RetOrNone<decltype(doTask(t,s))>
        {
        using Ret = RetOrNone<decltype(doTask(t,s))>;
        //return FixRet<Ret,Task&,Storage&>()(t,s);
        return std::declval<Ret>();
        }
    NoneType
    testRetImpl(Second)
        {
        return NoneType{};
        }
    };


//(4) Least preferred version which always compiles
template <typename Ret, typename Task, typename Storage>
Ret
callDoTask_Impl(Second, Task& t, Storage& s)
    {
    Error("1 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
//(3) This is preferred more than version (4) above
// since 0->int requires no cast.
// Template substitution will fail (not a compile-time error)
// if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret,typename Task, typename Storage>
auto 
callDoTask_Impl(First, Task& t, Storage& s)
    -> RetIfExists<Ret,decltype(doTask(t,s))>
    {
    return FixRet<Ret,decltype(doTask(t,s)),Task&,Storage&>()(t,s);
    }
//(2) Less preferred version which always compiles and passes
//the call on to callDoTask_Impl without the ManageStore& argument
template <typename Ret, typename Task, typename Storage>
Ret
callDoTask_Impl(Second, Task& t, Storage& s, ManageStore& m)
    {
    return callDoTask_Impl<Ret>(TryFirst,t,s);
    }
//(1) This is the preferred version since 0->int requires no cast
//Template substitution will fail (not a compile-time error)
//if doTask(t,d) is not defined ("SFINAE" trick)
template <typename Ret, typename Task, typename Storage>
auto 
callDoTask_Impl(First, Task& t, Storage& s, ManageStore& m)
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
    return callDoTask_Impl<Ret,Task,Storage>(TryFirst,t,s,m);
    }

/////////////////////////////////////////////////////

//(2-fail)
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Case2ConstNoMS(Task& t, const D& cd, ManageStore& m,Second)
    {
    if(!m.parg1().unique()) m.parg1() = m.parg1()->clone();
    auto* pd = static_cast<ITWrap<D>*>(m.parg1().get());
    return callDoTask<Ret>(t,pd->d,m);
    }
//(2-success)
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Case2ConstNoMS(Task& t, const D& cd, ManageStore& m,First)
    -> RetIfExists<Ret,decltype(doTask(t,cd))>
    {
    return FixRet<Ret,decltype(doTask(t,cd)),Task&,const D&>()(t,cd);
    }
//(1-fail)
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask_Case1ConstMS(Task& t, const D& cd, ManageStore& m,Second)
    {
    return cloneDoTask_Case2ConstNoMS<Ret>(t,cd,m,TryFirst);
    }
//(1-success) Preferred version since 0->int requires no conversion
//Attempts to call doTask(Task,const D&, ManageStore&) if defined
template<typename Ret, typename Task, typename D>
auto
cloneDoTask_Case1ConstMS(Task& t, const D& cd, ManageStore& m,First)
    -> RetIfExists<Ret,decltype(doTask(t,cd,m))>
    {
    return FixRet<Ret,decltype(doTask(t,cd,m)),Task&,const D&,ManageStore&>()(t,cd,m);
    }
template<typename Ret, typename Task, typename D>
Ret
cloneDoTask(Task& t, const D& cd, ManageStore& m)
    {
    return cloneDoTask_Case1ConstMS<Ret>(t,cd,m,TryFirst);
    }

/////////////////////

//(4) Least preferred version which always compiles
template <typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask_Impl(Second, Task& t, D1& d1, const D2& d2)
    {
    Error("2 parameter doTask not defined for specified task or data type");
    return Ret{};
    }
//(3) This is preferred more than version (4) above
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(First, Task& t, D1& d1, const D2& d2)
    -> RetIfExists<Ret,decltype(doTask(t,d1,d2))>
    {
    return FixRet<Ret,decltype(doTask(t,d1,d2)),Task&,D1&,const D2&>()(t,d1,d2);
    }
//(2) Less preferred version which always compiles and passes
//the call on to callDoTask_Impl (3) without the ManageStore& argument
template <typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask_Impl(Second, Task& t, D1& d1, const D2& d2, ManageStore& m)
    {
    return callDoTask_Impl<Ret>(TryFirst,t,d1,d2);
    }
//(1)
template <typename Ret, typename Task, typename D1, typename D2>
auto 
callDoTask_Impl(First, Task& t, D1& d1, const D2& d2, ManageStore& m)
    -> RetIfExists<Ret,decltype(doTask(t,d1,d2,m))>
    {
    return FixRet<Ret,decltype(doTask(t,d1,d2,m)),Task&,D1&,const D2&,ManageStore&>()(t,d1,d2,m);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
callDoTask(Task& t, D1& d1, const D2& d2, ManageStore& m)
    {
    //First try calling function labeled (1) above
    return callDoTask_Impl<Ret,Task,D1,D2>(TryFirst,t,d1,d2,m);
    }

/////////////////////

//(2-fail)
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_Case2ConstNoMS(Task& t, const D1& cd1, const D2& d2, ManageStore& m,Second) 
    {
    if(!m.parg1().unique()) m.parg1() = m.parg1()->clone();
    auto* pd1 = static_cast<std::remove_const_t<ITWrap<D1>>*>(m.parg1().get());
    return callDoTask<Ret>(t,pd1->d,d2,m);
    }
//(2-success)
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_Case2ConstNoMS(Task& t, const D1& cd1, const D2& d2, ManageStore& m,First) 
    -> RetIfExists<Ret,decltype(doTask(t,cd1,d2))>
    {
    return FixRet<Ret,decltype(doTask(t,cd1,d2)),Task&,const D1&,const D2&>()(t,cd1,d2);
    }
//(1-fail)
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask_Case1ConstMS(Task& t, const D1& cd1, const D2& d2, ManageStore& m,Second)
    {
    return cloneDoTask_Case2ConstNoMS<Ret,Task,D1,D2>(t,cd1,d2,m,0);
    }
//(1-success)
template<typename Ret, typename Task, typename D1, typename D2>
auto
cloneDoTask_Case1ConstMS(Task& t, const D1& cd1, const D2& d2, ManageStore& m,First) 
    -> RetIfExists<Ret,decltype(doTask(t,cd1,d2,m))>
    {
    return FixRet<Ret,decltype(doTask(t,cd1,d2,m)),Task&,const D1&,const D2&,ManageStore&>()(t,cd1,d2,m);
    }
template<typename Ret, typename Task, typename D1, typename D2>
Ret
cloneDoTask(Task& t, const D1& d1, const D2& d2, ManageStore& m)
    {
    return cloneDoTask_Case1ConstMS<Ret,Task,D1,D2>(t,d1,d2,m,TryFirst);
    }

void inline
check(const CPData& p)
    {
    if(!p) Error("doTask called on unallocated store pointer");
    }

} //namespace detail

template<typename Ret, typename Task, typename D>
Ret OneArg::
call(Task& t, const D& d, ManageStore& m)
    {
    return detail::callDoTask<Ret>(t,d,m);
    }

template<typename Ret, typename Task, typename D>
Ret OneArg::
call(Task& t, D& d, ManageStore& m)
    {
    return detail::cloneDoTask<Ret>(t,d,m);
    }

template<typename Ret, typename Task, typename D>
Ret TwoArgs::
call(Task& t, const D& d, ManageStore& m)
    {
    assert(m.hasArg2());
    CallWrap<Task,const D,Ret> w(t,d,m);
    m.arg2().plugInto(w);
    return w.getReturn();
    }

template<typename Ret, typename Task, typename D>
Ret TwoArgs::
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


//////
////// doTask methods
//////

template<typename T, typename... VArgs>
std::shared_ptr<ITData>
newITData(VArgs&&... vargs)
    {
    return std::make_shared<ITWrap<T>>(std::forward<VArgs>(vargs)...);
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
    RegisterTask<OneArg,Task,NoneType> r(std::forward<Task>(t));
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
    RegisterTask<OneArg,Task,NoneType> r(std::forward<Task>(t),&arg);
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
    RegisterTask<TwoArgs,Task,NoneType> r(std::forward<Task>(t),&arg1,&arg2);
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
    RegisterTask<TwoArgs,Task,NoneType> r(std::forward<Task>(t),&arg1,&arg2);
    arg1->plugInto(r);
    return std::move(r.getTask());
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

