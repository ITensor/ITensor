//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_TASK_TYPES_H_
#define __ITENSOR_TASK_TYPES_H_

#include "itensor/util/infarray.h"
#include "itensor/util/print.h"
#include "itensor/real.h"
#include "itensor/indexset.h"

namespace itensor {

//
// Task Types
// 

template<typename T>
struct Mult
    {
    T x;
    Mult(T x_) : x(x_) { }
    };

template<typename T>
const char*
typeNameOf(Mult<T> const&) { return "Mult"; }


struct MakeCplx { };

inline const char*
typeNameOf(MakeCplx const&) { return "MakeCplx"; }


struct GetElt
    {
    IndexSet const& is;
    IntArray const& inds;

    GetElt(IndexSet const& is_,
           IntArray const& inds_);
    };

inline const char*
typeNameOf(GetElt const&) { return "GetElt"; }

template<typename T>
struct SetElt
    {
    T elt;
    IndexSet const& is;
    IntArray const& inds;

    SetElt(T elt_,
           IndexSet const& is_,
           IntArray const& inds_);
    };

template<typename T>
const char*
typeNameOf(SetElt<T> const&) { return "SetElt"; }

struct NormNoScale { };

inline const char*
typeNameOf(NormNoScale const&) { return "NormNoScale"; }

struct PrintIT
    {
    std::ostream& s;
    LogNum const& x;
    IndexSet const& is;
    Real scalefac = 1.;
    bool print_data;

    PrintIT(std::ostream& s_,
            LogNum const& x_,
            IndexSet const& is_,
            bool print_data_)
        : s(s_), x(x_), is(is_), scalefac(1.), print_data(print_data_)
        { 
        if(!x.isTooBigForReal()) scalefac = x.real0();
        }

    template<typename D>
    void
    printInfo(D const& d, 
              std::string type_name,
              Real nrm_no_scale = -1)
        {
#ifndef USESCALE
        s << format("{norm=%.2f (%s)}\n",std::fabs(scalefac)*nrm_no_scale,type_name);
#else
        s << "{log(scale)=" << format("%.2f",x.logNum());
        if(nrm_no_scale > 0)
            {
            if(!x.isTooBigForReal()) s << ", norm=";
            else  s << ", norm(omitting large scale)=";
            s << format("%.2f",std::fabs(scalefac)*nrm_no_scale);
            }
        s << " (" << type_name << ")}\n";
#endif
        }
    };

inline const char*
typeNameOf(PrintIT const&) { return "PrintIT"; }

struct Conj { };

inline const char*
typeNameOf(Conj const&) { return "Conj"; }

struct CheckComplex { };

inline const char*
typeNameOf(CheckComplex const&) { return "CheckComplex"; }

struct SumEls
    {
    IndexSet const& is;
    SumEls(IndexSet const& is_) : is(is_) { }
    };

inline const char*
typeNameOf(SumEls const&) { return "SumEls"; }


template<typename F>
struct ApplyIT
    {
    F& f;

    ApplyIT(F&& f_) : f(f_)  { }

    template<typename T>
    void
    operator()(T& el)
        {
        applyITImpl(stdx::select_overload{},el,el);
        }

    template<typename T1, typename T2>
    void
    operator()(T1 from, T2& to)
        {
        applyITImpl(stdx::select_overload{},from,to);
        }

    private:

    template<typename T1, typename T2>
    void
    applyITImpl(stdx::choice<2>,T1, T2 &)
        {
        auto msg = format("Apply: function doesn't map %s->%s",typeName<T1>(),typeName<T2>());
        Error(msg);
        }
    template<typename T1, typename T2>
    auto
    applyITImpl(stdx::choice<1>,T1 from, T2 & to)
        -> stdx::if_compiles_return<void,decltype(to = f(from))>
        {
        to = f(from);
        }
    };

template<typename F>
const char*
typeNameOf(ApplyIT<F> const&) { return "ApplyIT"; }

namespace detail {

template<typename F, typename T>
void
applyType(F&&,T,long) { }

template<typename F, typename T>
auto
applyType(F&& f,T from,int)
    -> decltype(f(from))
    {
    return f(from);
    }
}

/// Traits for ApplyIT ///

template<typename F>
bool constexpr
realToCplx(ApplyIT<F> & A)
    {
    return std::is_same<decltype(detail::applyType(A.f,0.,0)),Cplx>::value;
    }

template<typename F>
bool constexpr
cplxToReal(ApplyIT<F> & A)
    {
    return std::is_same<decltype(detail::applyType(A.f,Cplx(0.),0)),Real>::value;
    }

template<typename T, typename F>
bool constexpr
switchesType(ApplyIT<F> & A)
    {
    return isReal<T>() ? realToCplx(A) : cplxToReal(A);
    }

template<typename T, typename F>
auto constexpr
resultTypeHelper(ApplyIT<F> const& A)
    -> stdx::conditional_t<std::is_void<decltype(detail::applyType(A.f,T{0.},0))>::value,
                           T,
                           decltype(detail::applyType(A.f,T{0.},0))>
    {
    using return_type = stdx::conditional_t<std::is_void<decltype(detail::applyType(A.f,T{0.},0))>::value,
                                            T,
                                            decltype(detail::applyType(A.f,T{0.},0))>;
    return return_type{0.};
    }

template<typename T,typename F>
using ApplyIT_result_of = decltype(resultTypeHelper<T>(std::declval<ApplyIT<F>>()));

///////////////////


template<typename F, typename T = typename std::invoke_result_t<F>>
struct GenerateIT
    {
    F& f;
    GenerateIT(F&& f_) : f(f_)  { }
    };

template<typename F, typename T>
const char*
typeNameOf(GenerateIT<F,T> const&) { return "GenerateIT"; }

template <typename F>
struct VisitIT
    {
    F& f;
    Real scale_fac = 1.;
    VisitIT(F&& f_, LogNum const& scale)
        : f(f_), scale_fac(scale.real0())
        { }
    };

template<typename F>
const char*
typeNameOf(VisitIT<F> const&) { return "VisitIT"; }

struct NNZBlocks { };

struct NNZ { };

struct IsDense { };

inline const char*
typeNameOf(NNZBlocks) { return "NNZBlocks"; }

template<typename T>
struct Fill
    {
    T x;
    Fill(T x_) : x(x_) { }
    };

template<typename T>
const char*
typeNameOf(Fill<T> const&) { return "Fill"; }

struct TakeReal { };
struct TakeImag { };

inline const char*
typeNameOf(TakeReal const&) { return "TakeReal"; }
inline const char*
typeNameOf(TakeImag const&) { return "TakeImag"; }

struct PlusEQ
    {
    using permutation = Permutation;
    private:
    const Permutation *perm_ = nullptr;
    const IndexSet *is1_ = nullptr,
                   *is2_ = nullptr;
    Real alpha_ = NAN;
    public:


    PlusEQ(Permutation const& P,
           IndexSet const& is1,
           IndexSet const& is2,
           Real alpha) :
        perm_(&P),
        is1_(&is1),
        is2_(&is2),
        alpha_(alpha)
        { }

    Real
    alpha() const { return alpha_; }

    Permutation const&
    perm() const { return *perm_; }

    IndexSet const&
    is1() const { return *is1_; }

    IndexSet const&
    is2() const { return *is2_; }
    };

inline const char*
typeNameOf(PlusEQ const&) { return "PlusEQ"; }

class Order
    {
    using permutation = Permutation;
    private:
    const Permutation *perm_ = nullptr;
    const IndexSet *is1_ = nullptr,
                   *is2_ = nullptr;
    public:


    Order(Permutation const& P,
          IndexSet const& is1,
          IndexSet const& is2) :
        perm_(&P),
        is1_(&is1),
        is2_(&is2)
        { }

    Permutation const&
    perm() const { return *perm_; }

    IndexSet const&
    is1() const { return *is1_; }

    IndexSet const&
    is2() const { return *is2_; }

    };

inline const char*
typeNameOf(Order const&) { return "Order"; }

#ifdef USESCALE
//
// Helper for Contract and NCProd
//
template<typename Storage>
Real
computeScalefac(Storage & dat)
    {
    //
    // TODO: better design could be
    //       to require each data type to
    //       an operator*= or similar
    //       Then, within e.g. ITReal make
    //       a TensorRef T and call norm(T)
    //       and T *= 1./scalefac
    //       Have TensorRef use dnrm2 and dscal
    auto scalefac = doTask(NormNoScale{},dat);
    //Here the NAN acts as a flag meaning "don't rescale"
    if(std::fabs(scalefac) < 1E-11) return NAN;
    doTask(Mult<Real>{1./scalefac},dat);
    return scalefac;
    }
#endif

struct Contract
    {
    IndexSet const& Lis;
    IndexSet const& Ris;
    IndexSet Nis; //new IndexSet
    Real scalefac = NAN;
    bool needresult = false;

    Contract(const IndexSet& Lis_,
             const IndexSet& Ris_)
      : Lis(Lis_),
        Ris(Ris_)
        { }

    Contract(const IndexSet& Lis_,
             const IndexSet& Ris_,
             const IndexSet& Nis_,
             bool needresult_ = false)
      : Lis(Lis_),
        Ris(Ris_),
        Nis(Nis_),
        needresult(needresult_)
        { }

    Contract(const Contract& other) = delete;
    Contract& operator=(const Contract& other) = delete;
    Contract(Contract&& other)
      : Lis(other.Lis),
        Ris(other.Ris),
        Nis(std::move(other.Nis)),
        scalefac(other.scalefac),
        needresult(other.needresult)
        { }

    };

inline const char*
typeNameOf(Contract const&) { return "Contract"; }

//Non-contracting product
struct NCProd
    {
    IndexSet const& Lis;
    IndexSet const& Ris;
    IndexSet Nis; //new IndexSet
    Real scalefac = NAN;

    NCProd(IndexSet const& Lis_,
           IndexSet const& Ris_)
      : Lis(Lis_),
        Ris(Ris_)
        { }

    NCProd(NCProd const& other) = delete;
    NCProd& operator=(NCProd const& other) = delete;
    NCProd(NCProd && other)
      : Lis(other.Lis),
        Ris(other.Ris),
        Nis(std::move(other.Nis)),
        scalefac(other.scalefac)
        { }
    };

inline const char*
typeNameOf(NCProd const&) { return "NCProd"; }

struct StorageType
    {
    enum Type
        { 
        Null=0, 
        DenseReal=1, 
        DenseCplx=2, 
        Combiner=3, 
        DiagReal=4, 
        DiagCplx=5,
        QDenseReal=6,
        QDenseCplx=7,
        QCombiner=8,
        QDiagReal=9,
        QDiagCplx=10,
        ScalarReal=11,
        ScalarCplx=12
        }; 
    };

inline const char*
typeNameOf(StorageType const&) { return "StorageType"; }


//template<typename T>
//void
//doTask(Write & W, T const& D)
//    {
//    auto sd = W.storageTypeDefined(stdx::select_overload{},D);
//    auto wd = W.writeDefined(stdx::select_overload{},D);
//    printfln("%d: storage type defined = %s, write defined = %s",typeNameOf(D),sd,wd);
//    W.writeType(stdx::select_overload{},D);
//    }


namespace detail {

void inline
checkEltInd(IndexSet const& is,
            IntArray const& inds)
    {
    for(auto k : range(inds))
        {
        auto i = inds[k];
        if(i < 0)
            {
            print("inds = ");
            for(auto j : inds) print(1+j," ");
            println();
            Error("Out of range: IndexVals are 1-indexed for getting tensor elements");
            }
        if(i >= dim(is[k]))
            {
            print("inds = ");
            for(auto j : inds) print(1+j," ");
            println();
            Error(format("Out of range: IndexVal at position %d has val %d > %s",1+k,1+i,Index(is[k])));
            }
        }
    }

} //namespace detail

inline GetElt::
GetElt(IndexSet const& is_,
       IntArray const& inds_)
  : is(is_),
    inds(inds_)
    { 
#ifdef DEBUG
    detail::checkEltInd(is,inds);
#endif
    }

template<typename T>
SetElt<T>::
SetElt(T elt_,
       IndexSet const& is_,
       IntArray const& inds_)
    : elt(elt_), is(is_), inds(inds_)
    { 
#ifdef DEBUG
    detail::checkEltInd(is,inds);
#endif
    }

struct CalcDiv 
    { 
    IndexSet const& is;
    CalcDiv(IndexSet const& is_) : is(is_) { }
    };

inline const char*
typeNameOf(CalcDiv const&) { return "CalcDiv"; }

struct IsEmpty { };
inline const char*
typeNameOf(IsEmpty const&) { return "IsEmpty"; }

template<typename V>
struct GetBlock
    {
    IndexSet const& is;
    Block const& block_ind;
    GetBlock(IndexSet const& is_,
             Block const& bi)
        : is(is_), block_ind(bi) { }
    };
inline const char*
typeNameOf(GetBlock<Real>) { return "GetBlock<Real>";}
inline const char*
typeNameOf(GetBlock<Cplx>) { return "GetBlock<Cplx>";}

struct RemoveQNs 
    {
    IndexSet const& is;
    RemoveQNs(IndexSet const& is_) : is(is_) {}
    };

inline const char*
typeNameOf(RemoveQNs) { return "RemoveQNs";}

struct ToDense
    {
    IndexSet const& is;
    ToDense(IndexSet const& is_) : is(is_) {}
    };

inline const char*
typeNameOf(ToDense) { return "ToDense";}

} //namespace itensor 

#endif
