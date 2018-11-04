#include "itensor/itensor_interface.h"

namespace itensor {


ITensor& ITensor::
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

Index
commonIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return Index{};
    }

Index
uniqueIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && !hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return Index{};
    }

ITensor
swapPrime(ITensor T, 
          int plev1, 
          int plev2,
          IndexType type)
    { 
    int tempLevel = 99999;
#ifdef DEBUG
    for(auto& I : T.inds())
        {
        if(I.primeLevel() == tempLevel) 
            {
            println("tempLevel = ",tempLevel);
            Error("swapPrime fails if an index has primeLevel==tempLevel");
            }
        }
#endif
    T.mapprime(plev1,tempLevel,type);
    T.mapprime(plev2,plev1,type);
    T.mapprime(tempLevel,plev2,type);
    return T; 
    }

Real
norm(ITensor const& T)
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

void
randomize(ITensor & T, Args const& args)
    {
    if(!T.store()) detail::allocReal(T);
#ifdef DEBUG
    if(!T) Error("default initialized tensor in randomize");
#endif
    auto cplx = args.getBool("Complex",false);
    if(cplx) T.generate(detail::quickranCplx);
    else     T.generate(detail::quickran);
    }

void ITensor::
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

ITensor
multSiteOps(ITensor A, ITensor const& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }

} //namespace itensor
