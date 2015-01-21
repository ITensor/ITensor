//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_IQTDATA_FUNCTIONS_H
#define __ITENSOR_IQTDATA_FUNCTIONS_H
#include "itdata/iqtdata.h"

namespace itensor {

template<typename T, int size>
struct IQGetElt
    {
    using BlockInd = IQIndexVal::BlockInd;
    using Inds = std::array<BlockInd,size>;

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
        long rank = is_.r();
        //Determine the number of the block containing
        //the element we are trying to get
        long block = 0;
        for(long j = 0, bstr = 1; j < rank; ++j)
            {
            block += (inds_[j].block-1)*bstr;
            bstr *= is_[j].nindex();
            }
        //
        //TODO: simplify above pattern (and other places) by create a generic
        //      object that wraps/adapts "Indexable" containers,
        //      applying a given lambda to them when operator[]
        //      is called, and supporting the size() function
        //      Maybe call it ContainerAdapter<Func> and have
        //      make_adapter(size,func)?
        //      (Pass all other args via the lambda.)
        //
        //Search for this block number to see if it's
        //contained in the list of offsets
        long offset = -1;
        for(const auto& io : d.offsets)
            if(io.ind == block) 
                {
                offset = io.offset;
                break;
                }
        //If not contained, element is zero
        //by quantum-number symmetry
        if(offset < 0)
            {
            elt_ = 0;
            return ITResult();
            }
        //Otherwise, move offset up to the 
        //position of the element in storage
        //and retrieve
        for(long j = 0, istr = 1; j < rank; ++j)
            {
            offset += (inds_[j].ind-1)*istr;
            istr *= is_[j].index(inds_[j].block).m();
            }
        if(offset >= d.data.size())
            {
            Error("IQIndexVal(s) out of range for IQTensor");
            }
        elt_ = d.data[offset];
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

