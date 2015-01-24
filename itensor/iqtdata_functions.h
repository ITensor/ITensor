//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_FUNCTIONS_H
#define __ITENSOR_IQTDATA_FUNCTIONS_H
#include "itdata/iqtdata.h"

namespace itensor {

template<typename T, size_t size>
struct IQGetElt
    {
    using Inds = std::array<long,size>;

    const IQIndexSet& is_;
    const Inds& inds_;
    T elt_;

    IQGetElt(const IQIndexSet& is,
             const Inds& inds)
        : 
        is_(is),
        inds_(inds)
        { }

    operator T() const { return elt_; }

    template <typename V,
              typename std::enable_if<std::is_convertible<V,T>::value>::type* = nullptr>
    ITResult
    operator()(const IQTData<V>& d)
        {
        auto* pelt = d.getElt(is_,inds_);
        if(pelt) elt_ = *pelt;
        else     elt_ = 0;
        return ITResult();
        }

    template <class D>
    ITResult
    operator()(const D& d)
        {
        throw ITError("IQTensor does not have requested element type");
        return ITResult();
        }
    };

template<size_t size>
struct IQSetEltReal
    {
    using Inds = std::array<long,size>;

    const IQIndexSet& is_;
    const Inds& inds_;
    Real elt_;

    IQSetEltReal(Real elt,
                 const IQIndexSet& is,
                 const Inds& inds)
        : 
        is_(is),
        inds_(inds),
        elt_(elt)
        { }

    template<typename T>
    ITResult
    operator()(IQTData<T>& d)
        {
        auto* pelt = d.getElt(is_,inds_);
        if(pelt) 
            {
            *pelt = T{elt_};
            }
        else 
            {
            Error("Setting IQTensor element non-zero would violate its symmetry.");
            }
        return ITResult();
        }

    template <class D>
    ITResult
    operator()(const D& d)
        {
        throw ITError("IQTensor does not have requested element type");
        return ITResult();
        }
    };

template <typename F>
class ApplyIQT
    {
    F& f_;
    public:
    ApplyIQT(F&& f) : f_(f) { }

    template <typename T,
              typename std::enable_if<std::is_same<T,std::result_of_t<F(T)>>::value>::type* = nullptr>
    ITResult
    operator()(IQTData<T>& d) const { return doApply(d); }

    private:

    template<typename T>
    ITResult
    doApply(T& d) const
        {
        for(auto& elt : d.data)
            elt = f_(elt);
        return ITResult();
        }
    };

template <typename F>
struct GenerateIQT
    {
    F& f_;
    public:
    GenerateIQT(F&& f)
        : f_(f)
        { }

    template <typename T>
    ITResult
    operator()(IQTData<T>& d) const { return doGen(d); }

    private:

    template<typename T>
    ITResult
    doGen(T& d) const
        {
        std::generate(d.data.begin(),d.data.end(),f_);
        return ITResult();
        }
    };

template <typename F>
class VisitIQT
    {
    F& f_;
    Real scale_fac;
    public:
    VisitIQT(F&& f, Real scale)
        : f_(f), scale_fac(scale)
        { }

    template <typename T>
    ITResult
    operator()(const T& d) const
        {
        for(const auto& elt : d.data)
            f_(elt*scale_fac);
        return ITResult();
        }
    };

}; //namespace itensor

#endif

