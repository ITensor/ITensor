#include "itensor/itensor_interface.h"

namespace itensor {

namespace detail {

void
allocReal(ITensor& T)
    {
    T.store() = newITData<DenseReal>(area(T.inds()),0);
    }

void
allocReal(IQTensor& T)
    {
    Error("Can't allocate IQTensor with undefined divergence");
    }

void
allocReal(ITensor& T, IntArray const& inds)
    {
    T.store() = newITData<DenseReal>(area(T.inds()),0);
    }

void
allocReal(IQTensor& T, IntArray const& inds)
    {
    QN div;
    for(size_t i = 0; i < T.inds().size(); ++i)
        {
        auto iv = (T.inds()[i])(1+inds[i]);
        div += iv.qn()*iv.index.dir();
        }
    T.store() = newITData<QDenseReal>(T.inds(),div);
    }

void
allocCplx(ITensor& T)
    {
    T.store() = newITData<DenseCplx>(area(T.inds()),0);
    }

void
checkArrows(IndexSet const& I1,
            IndexSet const& I2,
            bool shouldMatch = false)
    {
    //This space intentionally left blank
    }

void
checkArrows(IQIndexSet const& is1,
            IQIndexSet const& is2,
            bool shouldMatch = false)
    {
    for(auto I1 : is1)
    for(auto I2 : is2)
        {
        if(I1 == I2)
            {
            auto cond = shouldMatch ^ (I1.dir() == I2.dir());
            if(cond)
                {
                println("----------------------------------------");
                println("IQIndexSet 1 = \n",is1);
                println("----------------------------------------");
                println("IQIndexSet 2 = \n",is1);
                println("----------------------------------------");
                printfln("Mismatched IQIndex %s",I1);
                Error("Mismatched IQIndex arrows");
                }
            }
        }
    }

void
checkSameDiv(ITensor const& T1,
             ITensor const& T2)
    { }

void
checkSameDiv(IQTensor const& T1,
             IQTensor const& T2)
    {
    if(div(T1) != div(T2)) Error("div(T1) must equal div(T2) when adding T1+T2");
    }

} //namespace detail


template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
operator*=(ITensorT const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds());

    auto C = doTask(Contract<index_type>{L.inds(),R.inds()},
                    L.store(),
                    R.store());

    L.scale_ *= R.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }
template ITensorT<Index>& ITensorT<Index>::operator*=(ITensorT<Index> const& R);
template ITensorT<IQIndex>& ITensorT<IQIndex>::operator*=(ITensorT<IQIndex> const& R);

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
operator-=(const ITensorT& R)
    {
    auto& L = *this;
    if(&L == &R) 
        { 
        L.scale_ = LogNum(0); 
        return L;
        }
    L.scale().negate();
    L += R;
    L.scale().negate();
    return L;
    }
template ITensorT<Index>& ITensorT<Index>::operator-=(const ITensorT& R);
template ITensorT<IQIndex>& ITensorT<IQIndex>::operator-=(const ITensorT& R);

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
operator*=(Cplx z)
    {
    if(z.imag() == 0) return operator*=(z.real());
    doTask(Mult<Cplx>{z},store_);
    return *this;
    }
template ITensorT<Index>& ITensorT<Index>::operator*=(Cplx z);
template ITensorT<IQIndex>& ITensorT<IQIndex>::operator*=(Cplx z);

//Non-contracting product
template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
operator/=(ITensorT const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds(),true);

    auto C = doTask(NCProd<index_type>{L.inds(),R.inds()},
                    L.store(),
                    R.store());

    L.scale_ *= R.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }
template ITensorT<Index>& ITensorT<Index>::operator/=(ITensorT<Index> const& R);
template ITensorT<IQIndex>& ITensorT<IQIndex>::operator/=(ITensorT<IQIndex> const& R);

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
operator+=(ITensorT const& R)
    {
    auto& L = *this;
    if(!L) { return (L=R); }
    if(!R) Error("Right-hand-side of ITensor += is default constructed");
    if(&L == &R) return operator*=(2.);

    using permutation = typename PlusEQ<index_type>::permutation;

    auto P = permutation(L.inds().size());
#ifdef DEBUG
    try {
        calcPerm(R.inds(),L.inds(),P);
        }
    catch(std::exception const& e)
        {
        println("L = ",L);
        println("R = ",R);
        Error("ITensorT::operator+=: different index structure");
        }
#else
    calcPerm(R.inds(),L.inds(),P);
#endif

    if(Global::checkArrows()) 
        {
        auto shouldMatch = true;
        detail::checkArrows(L.inds(),R.inds(),shouldMatch);
        }

    //If L store unallocated, just assign from R
    if(!L.store() || L.scale().isZero()) return L = R;

#ifdef DEBUG
    detail::checkSameDiv(L,R);
#endif

    Real scalefac = 1.0;
    if(L.scale().magnitudeLessThan(R.scale())) 
        {
        L.scaleTo(R.scale()); 
        }
    else
        {
        scalefac = (R.scale()/L.scale()).real();
        }

    auto PEq = PlusEQ<index_type>{P,L.inds(),R.inds(),scalefac};
    doTask(PEq,L.store(),R.store());
    
    return L;
    } 
template ITensorT<Index>& ITensorT<Index>::operator+=(ITensorT<Index> const& R);
template ITensorT<IQIndex>& ITensorT<IQIndex>::operator+=(ITensorT<IQIndex> const& R);

template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
fill(Cplx z)
    {
    if(!store_) 
        {
        if(is_) detail::allocReal(*this);
        else Error("Can't fill default-constructed tensor");
        }
    scale_ = scale_type(1.);
    if(z.imag() == 0)
        doTask(Fill<Real>{z.real()},store_);
    else
        doTask(Fill<Cplx>{z},store_);
    return *this;
    }
template ITensorT<Index>& ITensorT<Index>::fill(Cplx z);
template ITensorT<IQIndex>& ITensorT<IQIndex>::fill(Cplx z);

template<typename IndexT> 
IndexT
commonIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return IndexT{};
    }
template Index commonIndex(const ITensorT<Index>& A, const ITensorT<Index>& B, IndexType t);
template IQIndex commonIndex(const ITensorT<IQIndex>& A, const ITensorT<IQIndex>& B, IndexType t);

template<typename IndexT> 
IndexT
uniqueIndex(const ITensorT<IndexT>& A, 
            const ITensorT<IndexT>& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && !hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return IndexT{};
    }
template Index uniqueIndex(const ITensorT<Index>& A, const ITensorT<Index>& B, IndexType t);
template IQIndex uniqueIndex(const ITensorT<IQIndex>& A, const ITensorT<IQIndex>& B, IndexType t);

template<typename IndexT>
ITensorT<IndexT>
swapPrime(ITensorT<IndexT> T, 
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
template ITensorT<Index> swapPrime(ITensorT<Index>, int, int, IndexType);
template ITensorT<IQIndex> swapPrime(ITensorT<IQIndex>, int, int, IndexType);

template<typename I>
Real
norm(ITensorT<I> const& T)
    {
#ifdef DEBUG
    if(!T) Error("Default initialized tensor in norm(ITensorT)");
#endif
    auto fac = std::fabs(T.scale().real0());
    return fac * doTask(NormNoScale{},T.store());
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
void
write(std::ostream& s, ITensorT<I> const& T)
    {
    write(s,T.inds());
    write(s,T.scale());
    auto type = StorageType::Null;
    if(T.store()) 
        {
        type = doTask(StorageType{},T.store());
        }
    write(s,type);
    if(T.store()) 
        {
        doTask(Write{s},T.store());
        }
    }
template void write(std::ostream& s, ITensorT<Index> const& T);
template void write(std::ostream& s, ITensorT<IQIndex> const& T);

template<class I>
ITensorT<I>
multSiteOps(ITensorT<I> A, ITensorT<I> const& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }
template ITensorT<Index> multSiteOps(ITensorT<Index> A, ITensorT<Index> const& B);
template ITensorT<IQIndex> multSiteOps(ITensorT<IQIndex> A, ITensorT<IQIndex> const& B);

} //namespace itensor
