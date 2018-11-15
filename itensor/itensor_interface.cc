#include "itensor/itensor_interface.h"

namespace itensor {


template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
fill(Cplx z)
    {
    if(!store_) 
        {
        if(is_) detail::allocReal(*this);
        else Error("Can't fill default-constructed tensor");
        }
    IF_USESCALE(scale_ = scale_type(1.);)
    if(z.imag() == 0)
        doTask(Fill<Real>{z.real()},store_);
    else
        doTask(Fill<Cplx>{z},store_);
    return *this;
    }
template ITensorT<Index>& ITensorT<Index>::fill(Cplx z);
template ITensorT<IQIndex>& ITensorT<IQIndex>::fill(Cplx z);

//TODO: make a version where the tags all have to match?
template<typename IndexT> 
IndexT
commonIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            TagSet const& t)
    {
    for(auto& I : A.inds())
        {
        if( (t == TagSet(All) || hasTags(I,t)) && hasIndex(B.inds(),I) ) 
            {
            return I;
            }
        }
    return IndexT{};
    }
template Index commonIndex(const ITensorT<Index>& A, const ITensorT<Index>& B, TagSet const& t);
template IQIndex commonIndex(const ITensorT<IQIndex>& A, const ITensorT<IQIndex>& B, TagSet const& t);

//TODO: make a version where the tags all have to match?
template<typename IndexT> 
IndexT
uniqueIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            TagSet const& t)
    {
    for(auto& I : A.inds())
        {
        if( (t == TagSet(All) || hasTags(I,t)) && !hasIndex(B.inds(),I) ) 
            {
            return I;
            }
        }
    return IndexT{};
    }
template Index uniqueIndex(const ITensorT<Index>& A, const ITensorT<Index>& B, TagSet const& t);
template IQIndex uniqueIndex(const ITensorT<IQIndex>& A, const ITensorT<IQIndex>& B, TagSet const& t);

template<typename I>
Real
norm(ITensorT<I> const& T)
    {
#ifdef DEBUG
    if(!T) Error("Default initialized tensor in norm(ITensorT)");
#endif

#ifndef USESCALE
    return doTask(NormNoScale{},T.store());
#else
    auto fac = std::fabs(T.scale().real0());
    return fac * doTask(NormNoScale{},T.store());
#endif
    }
template Real norm(ITensorT<Index> const& T);
template Real norm(ITensorT<IQIndex> const& T);

template<typename I>
void
randomize(ITensorT<I> & T, Args const& args)
    {
    if(!T.store()) detail::allocReal(T);
#ifdef DEBUG
    if(!T) Error("default initialized tensor in randomize");
#endif
    auto cplx = args.getBool("Complex",false);
    if(cplx) T.generate(detail::quickranCplx);
    else     T.generate(detail::quickran);
    }
template void randomize(ITensorT<Index> & T, Args const& args);
template void randomize(ITensorT<IQIndex> & T, Args const& args);


template<typename I>
void ITensorT<I>::
write(std::ostream& s) const
    {
    itensor::write(s,inds());
    itensor::write(s,scale());
    auto type = StorageType::Null;
    if(store()) 
        {
        type = doTask(StorageType{},store());
        }
    itensor::write(s,type);
    if(store()) 
        {
        doTask(Write{s},store());
        }
    }
template void ITensorT<Index>::write(std::ostream& s) const;
template void ITensorT<IQIndex>::write(std::ostream& s) const;

//
// TODO: should this be here?
// Maybe this is better for an example TEBD code?
//
template<class I>
ITensorT<I>
multSiteOps(ITensorT<I> A, ITensorT<I> const& B) 
    {
    A.prime("Site");
    A *= B;
    A.mapPrime(2,1,"Site");
    return A;
    }
template ITensorT<Index> multSiteOps(ITensorT<Index> A, ITensorT<Index> const& B);
template ITensorT<IQIndex> multSiteOps(ITensorT<IQIndex> A, ITensorT<IQIndex> const& B);

} //namespace itensor
