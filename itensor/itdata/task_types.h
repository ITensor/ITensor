//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
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

template<typename IndexT>
struct GetElt
    {
    using Inds = InfArray<long,28ul>;
    //using Inds = VarArray<long,31ul>;

    const IndexSetT<IndexT>& is;
    const Inds& inds;

    GetElt(const IndexSetT<IndexT>& is_,
           const Inds& inds_);
    };

template<typename T, typename IndexT>
struct SetElt
    {
    using Inds = InfArray<long,28ul>;
    //using Inds = VarArray<long,31ul>;
    T elt;
    const IndexSetT<IndexT>& is;
    const Inds& inds;

    SetElt(T elt_,
           const IndexSetT<IndexT>& is_,
           const Inds& inds_);
    };


struct NormNoScale { };

template<typename IndexT>
struct PrintIT
    {
    std::ostream& s;
    LogNum const& x;
    IndexSetT<IndexT> const& is;
    Real scalefac;
    bool print_data;

    PrintIT(std::ostream& s_,
            LogNum const& x_,
            IndexSetT<IndexT> const& is_,
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
        s << "{log(scale)=" << format("%.2f",x.logNum());
        if(nrm_no_scale > 0)
            {
            if(!x.isTooBigForReal()) s << ", norm=";
            else  s << ", norm(omitting large scale)=";
            s << format("%.2f",std::fabs(scalefac)*nrm_no_scale);
            }
        s << " (" << type_name << ")}\n";
        }

    void
    printVal(double val)
        {
        if(std::fabs(val) > 1E-10)
            s << val << "\n";
        else
            s << format("%.8E\n",val);
        }

    void
    printVal(Cplx const& val)
        {
        auto sgn = (val.imag() < 0 ? '-' : '+');
        if(std::norm(val) > 1E-10)
            {
            s << val.real() << sgn << std::fabs(val.imag()) << "i\n";
            }
        else
            {
            s << format("%.8E%s%.8E\n",val.real(),sgn,std::fabs(val.imag()));
            }
        }
    };

struct Conj { };

struct CheckComplex { };

template<typename IndexT>
struct SumEls
    {
    const IndexSetT<IndexT>& is;
    SumEls(IndexSetT<IndexT> const& is_) : is(is_) { }
    };

template<typename F>
struct ApplyIT
    {
    F& f;
    ApplyIT(F&& f_) : f(f_)  { }
    bool constexpr static 
    realToComplex()
        {
        return std::is_same<typename std::result_of<F(Real)>::type,Cplx>::value;
        }
    };

template<typename F, typename T = typename std::result_of<F()>::type>
struct GenerateIT
    {
    F& f;
    GenerateIT(F&& f_) : f(f_)  { }
    };

template <typename F>
struct VisitIT
    {
    F& f;
    Real scale_fac;
    VisitIT(F&& f_, const LogNum& scale)
        : f(f_), scale_fac(scale.real0())
        { }
    };

template<typename T>
struct Fill
    {
    T x;
    Fill(T x_) : x(x_) { }
    };

struct TakeReal { };
struct TakeImag { };

template<typename IndexT>
struct PlusEQ
    {
    using permutation = Permutation;
    using index_type = IndexT;
    using iset_type = IndexSetT<index_type>;
    private:
    const Permutation *perm_ = nullptr;
    const iset_type *is1_ = nullptr,
                    *is2_ = nullptr;
    Real fac_ = NAN;
    public:


    PlusEQ(Permutation const& P,
           iset_type const& is1,
           iset_type const& is2,
           Real fac) :
        perm_(&P),
        is1_(&is1),
        is2_(&is2),
        fac_(fac)
        { }

    Real
    fac() const { return fac_; }

    Permutation const&
    perm() const { return *perm_; }

    iset_type const&
    is1() const { return *is1_; }

    iset_type const&
    is2() const { return *is2_; }
    };

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
    if(std::fabs(scalefac) < 1E-11) return NAN;
    doTask(Mult<Real>{1./scalefac},dat);
    return scalefac;
    }

template<typename IndexT>
struct Contract
    {
    using index_type = IndexT;
    using iset_type = IndexSetT<IndexT>;

    iset_type const& Lis;
    iset_type const& Ris;
    iset_type Nis; //new IndexSet
    Real scalefac = NAN;
    bool needresult = false;

    Contract(const iset_type& Lis_,
             const iset_type& Ris_)
      : Lis(Lis_),
        Ris(Ris_)
        { }

    Contract(const iset_type& Lis_,
             const iset_type& Ris_,
             const iset_type& Nis_,
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

//Non-contracting product
template<typename IndexT>
struct NCProd
    {
    using index_type = IndexT;
    using iset_type = IndexSetT<IndexT>;

    iset_type const& Lis;
    iset_type const& Ris;
    iset_type Nis; //new IndexSet
    Real scalefac = NAN;

    NCProd(iset_type const& Lis_,
           iset_type const& Ris_)
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

struct StorageType
    {
    enum Type
        { 
        Null=0, 
        DenseReal=1, 
        DenseCplx=2, 
        ITCombiner=3, 
        DiagReal=4, 
        DiagCplx=5,
        IQTReal=6,
        IQTCombiner=7,
        IQTDiag=8
        }; 
    //const char*
    //name(Type t)
    //    {
    //    switch(t)
    //        {
    //        case DenseReal: return "DenseReal";
    //        case DenseCplx: return "DenseCplx";
    //        case ITCombiner: return "ICombiner";
    //        case DiagReal: return "DiagReal";
    //        case DiagCplx: return "DiagCplx";
    //        case IQTReal: return "IQTReal";
    //        case IQTReal: return "IQTReal";
    //        case IQTCombiner: return "IQTCombiner";
    //        case IQTDiag: return "IQTDiag";
    //        default: return "Null";
    //        }
    //    }
    //auto operator()(Dense<Real> const&) ->Type { return DenseReal_; }
    //auto operator()(Dense<Cplx> const&) ->Type { return DenseCplx_; }
    //auto operator()(ITCombiner  const&) ->Type { return ITCombiner_; }
    //auto operator()(Diag<Real>  const&) ->Type { return DiagReal_; }
    //auto operator()(Diag<Cplx>  const&) ->Type { return DiagCplx_; }
    //auto operator()(IQTReal     const&) ->Type { return IQTReal_; }
    //auto operator()(IQTCombiner const&) ->Type { return IQTCombiner_; }
    //auto operator()(IQTDiag     const&) ->Type { return IQTDiag_; }
    };

class Write
    {
    std::ostream& s;
    public:

    Write(std::ostream& s_) : s(s_) { }

    template<class T>
    void
    writeType(T const& data)
        {
        write(s,doTask(StorageType{},data));
        write(s,data); 
        }
    };

template<typename T>
void
doTask(Write & W, T const& D)
    {
    W.writeType(D);
    }


namespace detail {

template<typename I>
void
checkEltInd(const IndexSetT<I>& is,
            const typename GetElt<I>::Inds& inds)
    {
    for(size_t k = 0; k < inds.size(); ++k)
        {
        auto i = inds[k];
        if(i < 0)
            {
            print("inds = ");
            for(auto j : inds) print(1+j," ");
            println();
            Error("Out of range: IndexVals/IQIndexVals are 1-indexed for getting tensor elements");
            }
        if(i >= is[k].m())
            {
            print("inds = ");
            for(auto j : inds) print(1+j," ");
            println();
            Error(format("Out of range: IndexVal/IQIndexVal at position %d has val %d > %s",1+k,1+i,Index(is[k])));
            }
        }
    }

} //namespace detail

template<typename IndexT>
GetElt<IndexT>::
GetElt(const IndexSetT<IndexT>& is_,
       const Inds& inds_)
  : is(is_),
    inds(inds_)
    { 
#ifdef DEBUG
    detail::checkEltInd(is,inds);
#endif
    }

template<typename T, typename IndexT>
SetElt<T,IndexT>::
SetElt(T elt_,
       const IndexSetT<IndexT>& is_,
       const Inds& inds_)
    : elt(elt_), is(is_), inds(inds_)
    { 
#ifdef DEBUG
    detail::checkEltInd(is,inds);
#endif
    }

struct CalcDiv 
    { 
    IQIndexSet const& is;
    CalcDiv(IQIndexSet const& is_) : is(is_) { }
    };

} //namespace itensor 

#endif
