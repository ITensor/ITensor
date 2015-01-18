//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITDATA_FUNCTIONS_H
#define __ITENSOR_ITDATA_FUNCTIONS_H
#include "global.h"
#include "itdata/itdense.h"
#include "itdata/itdiag.h"
#include "itdata/itcombiner.h"
#include "indexset.h"
#include "simpletensor.h"
#include "contract.h"

namespace itensor {

template <typename F>
class ApplyIT
    {
    F& f_;
    public:
    ApplyIT(F&& f)
        : f_(f)
        { }

    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    ITResult
    operator()(ITDense<T>& d) const
        {
        for(auto& elt : d.data)
            {
            elt = f_(elt);
            }
        return ITResult();
        }

    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    ITResult
    operator()(ITDiag<T>& d) const
        {
        for(auto& elt : d.data)
            {
            elt = f_(elt);
            }
        return ITResult();
        }
    };

class Contract
    {
    const Label &Lind_,
                &Rind_;

    const IndexSet &Lis_,
                   &Ris_;

    //New IndexSet
    IndexSet Nis_;
    Real scalefac_ = -1;

    public:

    Contract(const IndexSet& Lis,
             const Label& Lind,
             const IndexSet& Ris,
             const Label& Rind)
        :
        Lind_(Lind),
        Rind_(Rind),
        Lis_(Lis),
        Ris_(Ris)
        { }

    IndexSet
    newIndexSet() { return std::move(Nis_); }
    Real
    scalefac() { return scalefac_; }

    ITResult
    operator()(const ITDense<Real>& a1,
               const ITDense<Real>& a2);

    ITResult
    operator()(const ITDense<Real>& d,
               const ITCombiner& C)
        {
        auto res = combine(d,Lis_,Ris_);
        if(!res) return ITResult::None;
        else     return std::move(res);
        }
    ITResult
    operator()(const ITCombiner& C,
               const ITDense<Real>& d)
        { 
        auto res = combine(d,Ris_,Lis_);
        if(!res) return ITResult::AssignPointer;
        else     return std::move(res);
        }

    ITResult
    operator()(const ITDiag<Real>& d,
               const ITDense<Real>& t)
        {
        return diagDense(d,Lis_,Lind_,t,Ris_,Rind_);
        }
    ITResult
    operator()(const ITDense<Real>& t,
               const ITDiag<Real>& d)
        { 
        return diagDense(d,Ris_,Rind_,t,Lis_,Lind_);
        }

    //ITResult
    //operator()(const ITDense<Real>& a1,
    //           const ITDense<Complex>& a2) const
    //    {
    //    ITDense<Complex> c1(a1);
    //    return operator()(c1,a2);
    //    }

    //ITResult
    //operator()(const ITDense<Complex>& a1,
    //           const ITDense<Real>& a2) const
    //    {
    //    ITDense<Complex> c2(a2);
    //    return operator()(a1,c2);
    //    }


    //template <typename T1, typename T2>
    //ITResult
    //operator()(const ITDense<T1>& a1,
    //           const ITDense<T2>& a2) const
    //    {
    //    using product_type = decltype(::std::declval<T1>() * ::std::declval<T2>());
    //    //static const auto One = product_type(1.),
    //    //                  Zero = product_type(0.);
    //    auto res = new ITDense<product_type>();
    //    //TODO:
    //    Error("Contract not implemented for tensors of different element types.");
    //    //btas::contract(One,a1.t_,Lind_,a2.t_,Rind_,Zero,res->t_,Nind_);
    //    return ITResult(res);
    //    }

    template <typename T1, typename T2>
    ITResult
    operator()(const T1& a1,const T2& a2) const
        {
        Error("Contract not implemented for this case");
        return ITResult();
        }

    private:

    enum SortOption { Sort, NoSort };
    void
    computeNis(SortOption sort);

    NewData
    combine(const ITDense<Real>& d,
            const IndexSet& dis,
            const IndexSet& Cis);

    ITResult
    diagDense(const ITDiag<Real>& d,
              const IndexSet& dis,
              const Label& dind,
              const ITDense<Real>& t,
              const IndexSet& tis,
              const Label& tind);
 
    };

class NormNoScale
    {
    Real nrm_;
    public:

    NormNoScale() : nrm_(0) { }

    operator Real() const { return nrm_; }

    template<typename T>
    ITResult
    operator()(const ITDense<T>& d) { return calc(d); }
    template<typename T>
    ITResult
    operator()(const ITDiag<T>& d) { return calc(d); }

    template<typename T>
    ITResult
    calc(const T& d)
        {
        for(const auto& elt : d.data)
            {
            nrm_ += std::norm(elt);
            }
        nrm_ = std::sqrt(nrm_);
        return ITResult();
        }
    };

class FillReal
    {
    Real r_;
    public:
    FillReal(Real r)
        : r_(r)
        { }

    ITResult
    operator()(ITDense<Real>& d) const;
    ITResult
    operator()(const ITDense<Complex>& d) const;
    ITResult
    operator()(ITDiag<Real>& d) const;
    ITResult
    operator()(const ITDiag<Complex>& d) const;

    template<typename T>
    ITResult
    operator()(const T& t) const
        {
        Error("FillReal: ITData type not implemented");
        return ITResult();
        }
    };

class FillCplx
    {
    Complex z_;
    public:
    FillCplx(Complex z)
        : z_(z)
        { }

    ITResult
    operator()(const ITDense<Real>& d) const;
    ITResult
    operator()(ITDense<Complex>& d) const;

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("Function not implemented."); return ITResult(); }
    };

template <typename F>
struct GenerateIT
    {
    F& f_;
    public:
    GenerateIT(F&& f)
        : f_(f)
        { }

    template <typename T>
    ITResult
    operator()(ITDense<T>& d) const { return doGen(d); }

    template <typename T>
    ITResult
    operator()(ITDiag<T>& d) const { return doGen(d); }

    private:

    template<typename T>
    ITResult
    doGen(T& d) const
        {
        std::generate(d.data.begin(),d.data.end(),f_);
        return ITResult();
        }
    };

template<typename T, int size>
struct GetElt
    {
    using Inds = std::array<long,size>;

    const IndexSet& is_;
    const Inds& inds_;
    T elt_;

    GetElt(const IndexSet& is,
           const Inds& inds)
        : 
        is_(is),
        inds_(inds)
        { }

    operator T() const { return elt_; }

    template <typename V,
              typename std::enable_if<std::is_convertible<V,T>::value>::type* = nullptr>
    ITResult
    operator()(const ITDense<V>& d)
        {
        elt_ = T{d.data[ind(is_,inds_)]};
        return ITResult();
        }

    template <typename V,
              typename std::enable_if<std::is_convertible<V,T>::value>::type* = nullptr>
    ITResult
    operator()(const ITDiag<V>& d)
        {
        auto first_i = (inds_.empty() ? 0 : inds_.front());
        //Check if inds_ reference an
        //element on the diagonal, else zero
        for(auto i : inds_)
            if(i != first_i)
                {
                elt_ = 0;
                return ITResult();
                }
        if(d.allSame())
            elt_ = d.val;
        else
            elt_ = d.data.at(first_i);
        return ITResult();
        }

    template <class D>
    ITResult
    operator()(const D& d)
        {
        throw ITError("ITensor does not have requested element type");
        return ITResult();
        }
    };

template<typename T, int size>
struct GetPtrElt
    {
    using Inds = std::array<long,size>;

    T* ptr_;
    const Inds& inds_;

    GetPtrElt(const Inds& inds)
        : inds_(inds)
        { }

    explicit operator T*() const { return ptr_; }

    template <typename V,
              typename std::enable_if<std::is_same<V,typename std::remove_const<T>::type>::value>::type* = nullptr>
    ITResult
    operator()(const ITDense<V>& d)
        {
        ptr_ = &(d.data.vref(d.data.ind(inds_)));
        return ITResult();
        }

    template <class D>
    ITResult
    operator()(const D& d)
        {
        throw ITError("ITensor does not have requested element type");
        return ITResult();
        }
    };

class MultComplex
    {
    Complex z_;
    public:
    MultComplex(Complex z) : z_(z) { }

    ITResult
    operator()(const ITDense<Real>& d) const;
    ITResult
    operator()(ITDense<Complex>& d) const;

    template<typename T>
    ITResult
    operator()(T& d) const { Error("MultComplex not defined for ITData type"); return ITResult(); }
    };

class MultReal
    {
    Real r_;
    public:
    MultReal(Real r)
        : r_(r)
        { }

    ITResult
    operator()(ITDense<Real>& d) const;
    ITResult
    operator()(ITDense<Complex>& d) const;

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("MultReal not implemented for ITData type."); return ITResult(); }
    };

class PlusEQ
    {
    Real fac_;
    const Permutation *P_ = nullptr;
    const IndexSet *is1_ = nullptr,
                   *is2_ = nullptr;
    bool permute_ = false;
    public:
    using permutation = Permutation;

    PlusEQ(Real fac)
        :
        fac_(fac)
        { }

    PlusEQ(const Permutation& P,
           const IndexSet& is1,
           const IndexSet& is2,
           Real fac)
        :
        fac_(fac),
        P_(&P),
        is1_(&is1),
        is2_(&is2),
        permute_(true)
        { }

    ITResult
    operator()(ITDense<Real>& a1,
               const ITDense<Real>& a2);

    ITResult
    operator()(ITDiag<Real>& a1,
               const ITDiag<Real>& a2);

    ITResult
    operator()(ITDense<Real>& a1,
               const ITDense<Complex>& a2)
        {
        Error("Real + Complex not implemented");
        //auto np = make_newdata<ITDense<Complex>>(a1);
        //operator()(*np,a2);
        //return ITResult(np);
        return ITResult();
        }

    ITResult
    operator()(ITDense<Complex>& a1,
               const ITDense<Real>& a2)
        {
        Error("Complex + Real not implemented");
        //ITDense<Complex> a2c(a2);
        //operator()(a1,a2c);
        return ITResult();
        }

    template <typename T1, typename T2>
    ITResult
    operator()(T1& a1,
               const T2& a2)
        {
        Error("Diag += not implemented");
        return ITResult();
        }
    };


struct PrintIT
    {
    std::ostream& s_;
    const LogNumber& x_;
    const IndexSet& is_;

    PrintIT(std::ostream& s,
            const LogNumber& x,
            const IndexSet& is)
        : s_(s), x_(x), is_(is)
        { }

    template<typename T>
    ITResult
    operator()(const ITDense<T>& d) const;

    template<typename T>
    ITResult
    operator()(const ITDiag<T>& d) const;

    ITResult
    operator()(const ITCombiner& d) const { s_ << " Combiner}\n"; return ITResult(); }

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("Function not implemented."); return ITResult(); }
    };

struct Read
    {
    std::istream& s_;
    Read(std::istream& s) : s_(s) { }
    
    template<typename DataType>
    ITResult
    operator()(DataType& d) const
        { 
        d.read(s_);
        return ITResult(); 
        }
    };

struct Write
    {
    std::ostream& s_;
    Write(std::ostream& s) : s_(s) { }
    
    template<typename DataType>
    ITResult
    operator()(const DataType& d) const
        { 
        d.write(s_);
        return ITResult(); 
        }
    };

class ReadWriteID
    {
    int id_ = 0;
    public:

    ReadWriteID() { }
    
    explicit operator int() const { return id_; }

    ITResult
    operator()(const ITDense<Real>& d) { id_ = 1; return ITResult(); }
    ITResult
    operator()(const ITDense<Complex>& d) { id_ = 2; return ITResult(); }

    ITResult static
    allocate(int id)
        {
        if(id == 1)
            return make_result<ITDense<Real>>();
        else 
        if(id == 2)
            return make_result<ITDense<Complex>>();
        else
            Error(format("ID %d not recognized",id));
        return ITResult();
        }
    };

template<long size>
class SetEltComplex
    {
    Complex elt_;
    const IndexSet& is_;
    const std::array<long,size>& inds_;
    public:
    SetEltComplex(Complex elt,
                  const IndexSet& is,
                  const std::array<long,size>& inds)
        : elt_(elt),
          is_(is),
          inds_(inds)
        { }

    ITResult
    operator()(const ITDense<Real>& d) const
        {
        auto nd = make_newdata<ITDense<Complex>>(d.data.cbegin(),d.data.cend());
        nd->data[ind(is_,inds_)] = elt_;
        return std::move(nd);
        }

    ITResult
    operator()(ITDense<Complex>& d) const
        {
        d.data[ind(is_,inds_)] = elt_;
        return ITResult();
        }

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("Function not implemented."); return ITResult(); }
    };

template<long size>
class SetEltReal
    {
    Real elt_;
    const IndexSet& is_;
    const std::array<long,size>& inds_;
    public:
    SetEltReal(Real elt,
               const IndexSet& is,
               const std::array<long,size>& inds)
        : elt_(elt),
          is_(is),
          inds_(inds)
        { }

    template<typename T>
    ITResult
    operator()(ITDense<T>& d) const
        {
        d.data[ind(is_,inds_)] = elt_;
        return ITResult();
        }

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("Function not implemented."); return ITResult(); }
    };

template <typename F>
class VisitIT
    {
    F& f_;
    Real scale_fac;
    public:
    VisitIT(F&& f, const LogNumber& scale)
        : f_(f), scale_fac(scale.real0())
        { }

    template <typename T>
    ITResult
    operator()(const T& d) const
        {
        for(const auto& elt : d.data)
            {
            f_(elt*scale_fac);
            }
        return ITResult();
        }
    };


}; //namespace itensor

#endif

