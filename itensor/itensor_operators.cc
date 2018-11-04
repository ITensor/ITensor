#include "itensor/itensor_interface.h"

namespace itensor {

namespace detail {

void
allocReal(ITensor& T)
    {
    if(hasQNs(T.inds())) Error("Can't allocate quantum ITensor with undefined divergence");
    T.store() = newITData<DenseReal>(area(T.inds()),0);
    }

void
allocReal(ITensor& T, IntArray const& inds)
    {
    if(!hasQNs(T.inds()))
        {
        T.store() = newITData<DenseReal>(area(T.inds()),0);
        }
    else
        {
        Error("allocReal not implemented for QN ITensor case");
        //QN div;
        //for(size_t i = 0; i < T.inds().size(); ++i)
        //    {
        //    auto iv = (T.inds()[i])(1+inds[i]);
        //    div += iv.qn()*iv.index.dir();
        //    }
        //T.store() = newITData<QDenseReal>(T.inds(),div);
        }
    }

void
allocCplx(ITensor& T)
    {
    if(hasQNs(T.inds())) Error("Can't allocate quantum ITensor with undefined divergence");
    T.store() = newITData<DenseCplx>(area(T.inds()),0);
    }


void
checkArrows(IndexSet const& is1,
            IndexSet const& is2,
            bool shouldMatch = false)
    {
    if(hasQNs(is1) && hasQNs(is2))
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
    }

void
checkSameDiv(ITensor const& T1,
             ITensor const& T2)
    {
    if(hasQNs(T1.inds()) && hasQNs(T2.inds()))
        {
        if(div(T1) != div(T2)) 
            {
            Error(format("div(T1)=%s must equal div(T2)=%s when adding T1+T2",div(T1),div(T2)));
            }
        }
    }

} //namespace detail


ITensor& ITensor::
operator*=(ITensor const& R)
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

    auto C = doTask(Contract{L.inds(),R.inds()},
                    L.store(),
                    R.store());

#ifdef USESCALE
    L.scale_ *= R.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;
#endif

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }

ITensor& ITensor::
order(IndexSet const& iset)
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
    auto bind = RangeBuilderT<IndexSet>(r);
    for(auto i : range(r))
        {
        bind.setIndex(P.dest(i),Ais[i]);
        }
    auto Bis = bind.build();

    auto O = Order{P,Ais,Bis};
    if(A.store())
        {
        doTask(O, A.store());
        }

    A.is_.swap(Bis);

    return A;
    }

#ifndef USESCALE

ITensor& ITensor::
operator*=(Real r)
    {
    doTask(Mult<Real>{r},store_);
    return *this;
    }

ITensor& ITensor::
operator/=(Real r)
    {
    auto fac = 1./r;
    doTask(Mult<Real>{fac},store_);
    return *this;
    }

#endif

ITensor& ITensor::
operator*=(Cplx z)
    {
    if(z.imag() == 0) return operator*=(z.real());
    doTask(Mult<Cplx>{z},store_);
    return *this;
    }

//Non-contracting product
ITensor& ITensor::
operator/=(ITensor const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds(),true);

    auto C = doTask(NCProd{L.inds(),R.inds()},
                    L.store(),
                    R.store());

#ifdef USESCALE
    L.scale_ *= R.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;
#endif

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }

void
daxpy(ITensor & L,
      ITensor const& R,
      Real alpha)
    {
    if(L.r() != R.r()) Error("ITensor::operator+=: different number of indices");

    using permutation = typename PlusEQ::permutation;

    auto P = permutation(L.inds().size());

    try {
        calcPerm(R.inds(),L.inds(),P);
        }
    catch(std::exception const& e)
        {
        println("L = ",L);
        println("R = ",R);
        Error("ITensoITensor::operator+=: different index structure");
        }

    if(Global::checkArrows()) 
        {
        auto shouldMatch = true;
        detail::checkArrows(L.inds(),R.inds(),shouldMatch);
        }

    if(!L.store()) Error("L not initialized in daxpy");

#ifdef DEBUG
    detail::checkSameDiv(L,R);
#endif

#ifdef USESCALE
    if(L.scale().magnitudeLessThan(R.scale())) 
        {
        L.scaleTo(R.scale()); 
        }
    else
        {
        alpha *= (R.scale()/L.scale()).real();
        }
#endif

    auto PEq = PlusEQ{P,L.inds(),R.inds(),alpha};
    doTask(PEq,L.store(),R.store());
    } 

ITensor& ITensor::
operator+=(ITensor const& R)
    {
    auto& L = *this;
    if(!L || !L.store()) { return (L=R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor += is default constructed");
    if(&L == &R) return operator*=(2.);

    daxpy(L,R,1.);
    
    return L;
    } 

ITensor& ITensor::
operator-=(ITensor const& R)
    {
    auto& L = *this;
    if(!L || !L.store()) { return (L = -R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor -= is default constructed");
    if(&L == &R) 
        { 
        L *= 0.;
        return L;
        }
    daxpy(L,R,-1.);
    return L;
    }



} //namespace itensor
