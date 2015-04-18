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
#include "itdata/iqtdata.h"
#include "indexset.h"
#include "simpletensor.h"
#include "contract.h"

namespace itensor {

template <typename F>
struct ApplyIT : RegisterFunc<ApplyIT<F>>
    {
    private:
    F& f_;
    public:
    ApplyIT(F&& f) : f_(f) { }

    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    void
    operator()(ITDense<T>& d) const { doApply(d); }
        
    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    void
    operator()(ITDiag<T>& d) const { doApply(d); }

    private:

    template<typename T>
    void
    doApply(T& d) const
        {
        for(auto& elt : d.data)
            elt = f_(elt);
        }
    };





template <typename F>
struct GenerateIT : RegisterFunc<GenerateIT<F>>
    {
    private:
    F& f_;
    public:
    GenerateIT(F&& f) : f_(f) { }

    template <typename T>
    void
    operator()(ITDense<T>& d) const { doGen(d); }

    template <typename T>
    void
    operator()(ITDiag<T>& d) const { doGen(d); }

    private:

    template<typename T>
    void
    doGen(T& d) const
        {
        std::generate(d.data.begin(),d.data.end(),f_);
        }
    };

template<typename T, int size>
struct GetElt : RegisterFunc<GetElt<T,size>>
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
    void
    operator()(const ITDense<V>& d)
        {
        elt_ = T{d.data[ind(is_,inds_)]};
        }

    template <typename V,
              typename std::enable_if<std::is_convertible<V,T>::value>::type* = nullptr>
    void
    operator()(const ITDiag<V>& d)
        {
        auto first_i = (inds_.empty() ? 0 : inds_.front());
        //Check if inds_ reference an
        //element on the diagonal, else zero
        for(auto i : inds_)
            if(i != first_i)
                {
                elt_ = 0;
                return;
                }
        if(d.allSame())
            elt_ = d.val;
        else
            elt_ = d.data.at(first_i);
        }

    //template <class D>
    //ITResult
    //operator()(const D& d)
    //    {
    //    throw ITError("ITensor does not have requested element type");
    //    return ITResult();
    //    }
    };

template<typename T, int size>
struct GetPtrElt : RegisterFunc<GetPtrElt<T,size>>
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
    void
    operator()(const ITDense<V>& d)
        {
        ptr_ = &(d.data.vref(d.data.ind(inds_)));
        }

    template <class D>
    void
    operator()(const D& d)
        {
        throw ITError("ITensor does not have requested element type");
        }
    };

struct MultComplex : RegisterFunc<MultComplex>
    {
    private:
    Complex z_;
    public:
    MultComplex(Complex z) : z_(z) { }

    void
    operator()(const ITDense<Real>& d);
    void
    operator()(ITDense<Complex>& d);

    template<typename T>
    void
    operator()(T& d) const { Error("MultComplex not defined for ITData type"); }
    };



struct PrintIT : RegisterFunc<PrintIT>
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
    void
    operator()(const ITDense<T>& d) const;

    template<typename T>
    void
    operator()(const ITDiag<T>& d) const;

    void
    operator()(const ITCombiner& d) const { s_ << " Combiner}\n"; }
    };

//struct Read : RegisterFunc<Read>
//    {
//    std::istream& s_;
//    Read(std::istream& s) : s_(s) { }
//    
//    template<typename DataType>
//    void
//    operator()(DataType& d) const
//        { 
//        d.read(s_);
//        }
//    };
//
//struct Write : RegisterFunc<Write>
//    {
//    std::ostream& s_;
//    Write(std::ostream& s) : s_(s) { }
//    
//    template<typename DataType>
//    void
//    operator()(const DataType& d) const
//        { 
//        d.write(s_);
//        }
//    };

template<size_t size>
struct SetEltComplex : RegisterFunc<SetEltComplex<size>>
    {
    private:
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

    void
    operator()(const ITDense<Real>& d)
        {
        auto nd = this->template setNewData<ITDense<Complex>>(d.data.cbegin(),d.data.cend());
        nd->data[ind(is_,inds_)] = elt_;
        }

    void
    operator()(ITDense<Complex>& d)
        {
        d.data[ind(is_,inds_)] = elt_;
        }
    };

template<long size>
struct SetEltReal : RegisterFunc<SetEltReal<size>>
    {
    private:
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
    void
    operator()(ITDense<T>& d) const
        {
        d.data[ind(is_,inds_)] = elt_;
        }
    };

template <typename F>
struct VisitIT : RegisterFunc<VisitIT<F>>
    {
    private:
    F& f_;
    Real scale_fac;
    public:
    VisitIT(F&& f, const LogNumber& scale)
        : f_(f), scale_fac(scale.real0())
        { }

    template <typename T>
    void
    operator()(const T& d) const
        {
        for(const auto& elt : d.data)
            {
            f_(elt*scale_fac);
            }
        }
    };


}; //namespace itensor

#endif

