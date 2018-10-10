//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SCALAR_H
#define __ITENSOR_SCALAR_H

#include "itensor/itdata/task_types.h"
//#include "itensor/itdata/itdata.h"
#include "itensor/itdata/dotask.h"
#include "itensor/iqindex.h"
#include "itensor/detail/call_rewrite.h"

namespace itensor {

template<typename T>
class Scalar;

using ScalarCplx = Scalar<Cplx>;
using ScalarReal = Scalar<Real>;

template<typename T>
class Scalar
    {
    static_assert(not std::is_const<T>::value,
                  "Template argument to Scalar storage should not be const");
    public:
    using value_type = T;

    value_type val = 0;

    Scalar() { }

    explicit
    Scalar(value_type z) : val(z) { }

    value_type*
    data() { return &val; }
    value_type const*
    data() const { return &val; }
    };

namespace detail {
Real inline
ensureReal(Cplx z) { return z.real(); }
Real inline
ensureReal(Real r) { return r; }
void inline
assign(Cplx & l, Real r) { l = r; }
void inline
assign(Cplx & l, Cplx r) { l = r; }
void inline
assign(Real & l, Real r) { l = r; }
void inline
assign(Real & l, Cplx r) { Error("Cannot assign Real = Cplx"); }
} //namespace detail

const char*
typeNameOf(ScalarReal const& d);
const char*
typeNameOf(ScalarCplx const& d);

template<typename T>
bool constexpr
isReal(Scalar<T> const& t) { return std::is_same<T,Real>::value; }

template<typename T>
bool constexpr
isCplx(Scalar<T> const& t) { return std::is_same<T,Cplx>::value; }

template<typename T>
void 
read(std::istream& s, Scalar<T> & d)
    {
    itensor::read(s,d.val);
    }

template<typename T>
void
write(std::ostream& s, Scalar<T> const& d)
    {
    itensor::write(s,d.val);
    }


template<typename F, typename T>
void
doTask(ApplyIT<F>& A, Scalar<T> & d, ManageStore & m)
    { 
    using new_type = ApplyIT_result_of<T,F>;
    if(switchesType<T>(A))
        {
        auto *nd = m.makeNewData<Scalar<new_type>>();
        A(d.val,nd->val);
        }
    else
        {
        auto *md = m.modifyData(d);
        A(md->val);
        }
    }

template<typename F, typename T>
void
doTask(VisitIT<F>& V, Scalar<T> const& d)
    { 
    detail::call<void>(V.f,V.scale_fac * d.val);
    }

template<typename I, typename T>
Cplx 
doTask(GetElt<I> const& g, Scalar<T> const& d) { return d.val; }

template<typename E, typename I, typename T>
void
doTask(SetElt<E,I> const& S, Scalar<T> const& d, ManageStore & m)
    {
    if(not std::is_same<E,T>::value)
        {
        auto& nd = *m.makeNewData<Scalar<E>>();
        nd.val = S.elt;
        }
    else
        {
        auto& dnc = *m.modifyData(d);
        detail::assign(dnc.val,S.elt);
        }
    }

template<typename E, typename T>
void
doTask(Fill<E> const& f, ScalarReal const& d, ManageStore & m)
    {
    if(std::is_same<E,T>::value)
        {
        auto& dnc = *m.modifyData(d);
        dnc.val = f.x;
        }
    else
        {
        m.makeNewData<Scalar<E>>(f.x);
        }
    }

template<typename E, typename T>
void
doTask(Mult<E> const& M, Scalar<T> const& d, ManageStore & m)
    {
    if(std::is_same<E,Cplx>::value && isReal(d))
        {
        m.makeNewData<ScalarCplx>(d.val*M.x);
        }
    else
        {
        auto& dref = *m.modifyData(d);
        detail::assign(dref.val,dref.val*M.x);
        }
    }

template<typename T>
Real
doTask(NormNoScale, Scalar<T> const& d) { return std::abs(d.val); }

void inline
doTask(Conj, ScalarReal const& d) { }

void inline
doTask(Conj, ScalarCplx & d) { std::conj(d.val); }

void inline
doTask(TakeReal, ScalarReal const& ) { /*nothing to do*/ }

void inline
doTask(TakeReal, ScalarCplx const& d, ManageStore & m)
    {
    m.makeNewData<ScalarReal>(d.val.real());
    }

void inline
doTask(TakeImag, ScalarReal & d) { d.val = 0.; }

void inline
doTask(TakeImag, ScalarCplx const& d, ManageStore & m)
    {
    m.makeNewData<ScalarReal>(d.val.imag());
    }

void inline
doTask(MakeCplx, ScalarReal const& d, ManageStore & m)
    {
    m.makeNewData<ScalarCplx>(d.val);
    }

void inline
doTask(MakeCplx, ScalarCplx const& d) { /*nothing to do*/ }

template<typename T>
bool constexpr
doTask(CheckComplex, Scalar<T> const& d) { return isCplx(d); }

template<typename I, typename T>
void
doTask(PrintIT<I>& P, Scalar<T> const& d)
    {
    auto name = std::is_same<T,Real>::value ? "Scalar Real"
                                            : "Scalar Cplx";
    P.printInfo(d,name,doTask(NormNoScale{},d));
    if(P.print_data)
        {
        P.s << "  " << formatVal(P.scalefac*d.val) << "\n";
        }
    }

template<typename I, typename T>
Cplx
doTask(SumEls<I>, Scalar<T> const& d) { return d.val; }

auto constexpr inline
doTask(StorageType const& S, ScalarReal const& d) ->StorageType::Type { return StorageType::ScalarReal; }

auto constexpr inline
doTask(StorageType const& S, ScalarCplx const& d) ->StorageType::Type { return StorageType::ScalarCplx; }


//
// These Contract implementations are correct, but currently
// not used due to rank==0 case caught by ITensorT operator*= function
//
//template<typename I, typename T, typename StoreType,
//         class = typename stdx::enable_if_t<containsType<StorageTypes,stdx::decay_t<StoreType>>{}> >
//void
//doTask(Contract<I> & C,
//       Scalar<T> const& L,
//       StoreType const& R,
//       ManageStore & m)
//    {
//    C.Nis = C.Ris;
//    if(isReal(L))
//        {
//        //Make left-hand tensor have R as storage
//        m.assignPointerRtoL();
//        C.scalefac = detail::ensureReal(L.val);
//        }
//    else //Scalar L is complex
//        {
//        m.makeNewData<StoreType>(R);
//        doTask(Mult<Cplx>(L.val),m.newData());
//        }
//    }
//
//template<typename I, typename T,typename StoreType,
//         class = typename stdx::enable_if_t<containsType<StorageTypes,stdx::decay_t<StoreType>>{}> >
//void
//doTask(Contract<I> & C,
//       StoreType const& L,
//       Scalar<T> const& R,
//       ManageStore & m)
//    {
//    C.Nis = C.Lis;
//    if(isReal(R))
//        {
//        //Just multiply value of Scalar into scalefac 
//        //of Contract task object
//        C.scalefac = detail::ensureReal(R.val);
//        }
//    else //Scalar R is complex
//        {
//        m.modifyData();
//        doTask(Mult<Cplx>(R.val),m.parg1());
//        }
//    }

template<typename I, typename T1, typename T2>
void
doTask(PlusEQ<I> const& P,
       Scalar<T1> const& d1,
       Scalar<T2> const& d2,
       ManageStore & m)
    {
    auto s = d1.val + P.alpha()*d2.val;
    if(isReal(d1) && isCplx(d2))
        {
        m.makeNewData<ScalarCplx>(s);
        }
    else
        {
        auto& d1ref = *m.modifyData(d1);
        detail::assign(d1ref.val,s);
        }
    }

template<typename T>
QN
doTask(CalcDiv const& c,Scalar<T> const& d) { return QN(); }


} //namespace itensor

#endif
