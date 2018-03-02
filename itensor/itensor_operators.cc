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
                println("IQIndexSet 2 = \n",is2);
                println("----------------------------------------");
                printfln("Mismatched IQIndex from set 1 %s",I1);
                printfln("Mismatched IQIndex from set 2 %s",I2);
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
    if(div(T1) != div(T2)) 
        {
        Error(format("div(T1)=%s must equal div(T2)=%s when adding T1+T2",div(T1),div(T2)));
        }
    }

} //namespace detail


template<typename IndexT> 
ITensorT<IndexT>& ITensorT<IndexT>::
operator*=(ITensorT const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(L.r() == 0)
        {
        auto z = L.cplx();
        *this = R*z;
        return *this;
        }
    else if(R.r()==0)
        {
        auto z = R.cplx();
        *this *= z;
        return *this;
        }

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
order(IndexSetT<IndexT> const& iset)
    {
    auto& A = *this;
    auto Ais = A.inds();
    auto r = Ais.r();

    if(size_t(r) != size_t(iset.r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",Ais,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n",iset,"\n");
        println("---------------------------------------------");
        Error(format("Wrong number of Indexes passed to order (expected %d, got %d)",r,iset.r()));
        }

    // Get permutation
    auto P = Permutation(r);
    calcPerm(Ais,iset,P);
    if(isTrivial(P))
        {
        return A;
        }
    // If not trivial, use permutation to get new index set
    // This is necessary to preserve the proper arrow direction of IQIndex
    auto bind = RangeBuilderT<IndexSetT<IndexT>>(r);
    for(auto i : range(r))
        {
        bind.setIndex(P.dest(i),Ais[i]);
        }
    auto Bis = bind.build();

    auto O = Order<IndexT>{P,Ais,Bis};
    if(A.store())
        doTask(O, A.store());

    A.is_.swap(Bis);

    return A;
    }
template ITensorT<Index>& ITensorT<Index>::order(IndexSetT<Index> const& iset);
template ITensorT<IQIndex>& ITensorT<IQIndex>::order(IndexSetT<IQIndex> const& iset);

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
    if(!L) { return (L=R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor += is default constructed");
    if(&L == &R) return operator*=(2.);

    if(this->r() != R.r()) Error("ITensorT::operator+=: different number of indices");

    using permutation = typename PlusEQ<index_type>::permutation;

    auto P = permutation(L.inds().size());

    try {
        calcPerm(R.inds(),L.inds(),P);
        }
    catch(std::exception const& e)
        {
        println("L = ",L);
        println("R = ",R);
        Error("ITensorT::operator+=: different index structure");
        }

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


} //namespace itensor
