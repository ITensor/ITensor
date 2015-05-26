//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TASK_TYPES_H_
#define __ITENSOR_TASK_TYPES_H_

#include "itensor/util/infarray.h"
#include "itensor/tensor/permute.h"

namespace itensor {

//
// Task Types
// 

template<size_t size>
struct GetElt
    {
    using Inds = std::array<long,size>;

    const IndexSet& is;
    const Inds& inds;

    GetElt(const IndexSet& is_,
           const Inds& inds_)
      : is(is_),
        inds(inds_)
        { 
#ifdef DEBUG
        for(size_t k = 0; k < inds.size(); ++k)
            {
            auto i = inds[k];
            if(i < 0)
                {
                print("inds = ");
                for(auto j : inds) print(1+j," ");
                println();
                Error("Out of range: IndexVals are 1-indexed for getting ITensor elements");
                }
            if(i >= is[k].m())
                {
                print("inds = ");
                for(auto j : inds) print(1+j," ");
                println();
                Error(format("Out of range: IndexVal at position %d has val %d > %d",1+k,1+i,is[k].m()));
                }
            }
#endif
        }
    };

template<typename T, size_t size>
struct SetElt
    {
    using Inds = std::array<long,size>;
    T elt;
    const IndexSet& is;
    const Inds& inds;

    SetElt(T elt_,
           const IndexSet& is_,
           const Inds& inds_)
        : elt(elt_), is(is_), inds(inds_)
        { }
    };

struct MultReal
    {
    Real r;
    MultReal(Real r_) : r(r_) { }
    };

struct MultCplx
    {
    Cplx z;
    MultCplx(Cplx z_) : z(z_) { }
    };

struct NormNoScale
    {
    const IndexSet& is;
    NormNoScale(const IndexSet& is_) : is(is_) { }
    };

struct PrintIT
    {
    std::ostream& s;
    const LogNumber& x;
    const IndexSet& is;
    Real scalefac;
    bool print_data;

    PrintIT(std::ostream& s_,
            const LogNumber& x_,
            const IndexSet& is_,
            bool print_data_)
        : s(s_), x(x_), is(is_), scalefac(1.), print_data(print_data_)
        { 
        if(!x.isTooBigForReal()) scalefac = x.real0();
        }

    template<typename D>
    void
    printInfo(const D& d, 
              std::string type_name,
              Real nrm_no_scale = -1)
        {
        s << "{log(scale)=" << format("%.2f",x.logNum());
        if(nrm_no_scale > 0)
            {
            if(!x.isTooBigForReal()) s << ", norm=";
            else  s << ", norm(omitting large scale)=";
            s << format("%.2f",fabs(scalefac)*nrm_no_scale);
            }
        s << " (" << type_name << ")}\n";
        }
    };

struct Conj { };

struct CheckComplex { };

struct SumEls
    {
    const IndexSet& is;
    SumEls(const IndexSet& is_) : is(is_) { }
    };

template<typename F>
struct ApplyIT
    {
    F& f;
    ApplyIT(F&& f_) : f(f_)  { }
    };

template<typename F>
struct GenerateIT
    {
    F& f;
    GenerateIT(F&& f_) : f(f_)  { }
    bool constexpr static 
    isComplex()
        {
        return std::is_same<typename std::result_of<F()>::type,Cplx>::value;
        }
    };

template <typename F>
struct VisitIT
    {
    F& f;
    Real scale_fac;
    VisitIT(F&& f_, const LogNumber& scale)
        : f(f_), scale_fac(scale.real0())
        { }
    };

struct FillReal
    {
    Real r;
    FillReal(Real r_) : r(r_) { }
    };

struct FillCplx
    {
    Cplx z;
    FillCplx(Cplx z_) : z(z_) { }
    };

struct TakeReal { };
struct TakeImag { };

struct Contract
    {
    const Label &Lind,
                &Rind;

    const IndexSet &Lis,
                   &Ris;

    IndexSet Nis; //new IndexSet
    Real scalefac = NAN;

    Contract(const IndexSet& Lis_,
             const Label& Lind_,
             const IndexSet& Ris_,
             const Label& Rind_) :
        Lind(Lind_),
        Rind(Rind_),
        Lis(Lis_),
        Ris(Ris_)
        { }
    };


} //namespace itensor 

#endif
