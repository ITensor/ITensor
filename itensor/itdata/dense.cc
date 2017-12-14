//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/dense.h"
#include "itensor/itdata/itdata.h"
#include "itensor/itdata/itlazy.h"
#include "itensor/indexset.h"
#include "itensor/util/range.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/lapack_wrap.h"

namespace itensor {

const char*
typeNameOf(DenseReal const& d) { return "DenseReal"; }
const char*
typeNameOf(DenseCplx const& d) { return "DenseCplx"; }

Cplx 
doTask(GetElt<Index> const& g, DenseReal const& d)
    {
    return Cplx(d[offset(g.is,g.inds)],0.);
    }

Cplx 
doTask(GetElt<Index> const& g, DenseCplx const& d)
    {
    return d[offset(g.is,g.inds)];
    }

template<typename E, typename T>
struct SetEltHelper
    {
    void static
    set(SetElt<E,Index> const& S, Dense<T> const& D, ManageStore& m)
        {
        auto& Dnc = *m.modifyData(D);
        Dnc[offset(S.is,S.inds)] = S.elt;
        }
    };
template<>
struct SetEltHelper<Cplx,Real>
    {
    void static
    set(SetElt<Cplx,Index> const& S, DenseReal const& D, ManageStore & m)
        {
        auto& nd = *m.makeNewData<DenseCplx>(D.begin(),D.end());
        nd[offset(S.is,S.inds)] = S.elt;
        }
    };

template<typename E, typename T>
void
doTask(SetElt<E,Index> const& S, Dense<T> const& D, ManageStore & m)
    {
    SetEltHelper<E,T>::set(S,D,m);
    }
template
void
doTask(SetElt<Real,Index> const& S, DenseReal const& D, ManageStore & m);
template
void
doTask(SetElt<Real,Index> const& S, DenseCplx const& D, ManageStore & m);
template
void
doTask(SetElt<Cplx,Index> const& S, DenseReal const& D, ManageStore & m);
template
void
doTask(SetElt<Cplx,Index> const& S, DenseCplx const& D, ManageStore & m);


void
doTask(Fill<Real> const& f, DenseReal & D)
    {
    stdx::fill(D,f.x);
    }
void
doTask(Fill<Real> const& f, DenseCplx const& D, ManageStore & m)
    {
    m.makeNewData<DenseReal>(D.size(),f.x);
    }
void
doTask(Fill<Cplx> const& f, DenseReal const& D, ManageStore & m)
    {
    m.makeNewData<DenseCplx>(D.size(),f.x);
    }
void
doTask(Fill<Cplx> const& f, DenseCplx & D)
    {
    stdx::fill(D,f.x);
    }


void
doTask(Mult<Cplx> const& M, Dense<Cplx> & D)
    {
    for(auto& el : D) el *= M.x;
    }
void
doTask(Mult<Cplx> const& M, Dense<Real> const& D, ManageStore & m)
    {
    auto nd = m.makeNewData<DenseCplx>(D.begin(),D.end());
    doTask(M,*nd);
    }

template<typename T>
void
doTask(Mult<Real> const& M, Dense<T> & D)
    {
    auto d = realData(D);
    dscal_wrapper(d.size(),M.x,d.data());
    }
template
void
doTask(Mult<Real> const& M, DenseReal & D);
template
void
doTask(Mult<Real> const& M, DenseCplx & D);

template<typename T>
Real
doTask(NormNoScale, Dense<T> const& D) 
    { 
    auto d = realData(D);
    return dnrm2_wrapper(d.size(),d.data());
    }
template
Real
doTask(NormNoScale, DenseReal const& D);
template
Real
doTask(NormNoScale, DenseCplx const& D);

void
doTask(Conj,DenseReal const& D) { /*Nothing to conj*/ }

void
doTask(Conj,DenseCplx & D) 
    { 
    for(auto& el : D) el = std::conj(el);
    }

void
doTask(TakeReal, DenseReal const& D) { /*Already real*/ }

void
doTask(TakeImag, DenseReal & D)
    { 
    //Set all elements to zero
    doTask(Mult<Real>{0.},D);
    }

void
doTask(TakeReal, DenseCplx const& D, ManageStore & m) 
    { 
    auto& nD = *m.makeNewData<DenseReal>(D.size());
    for(auto n : range(D.size()))
        {
        nD[n] = D[n].real();
        }
    }

void
doTask(TakeImag, DenseCplx const& D, ManageStore & m) 
    { 
    auto& nD = *m.makeNewData<DenseReal>(D.size());
    for(auto n : range(D.size()))
        {
        nD[n] = D[n].imag();
        }
    }

void
doTask(MakeCplx, DenseReal const& d, ManageStore & m)
    {
    m.makeNewData<DenseCplx>(d.begin(),d.end());
    }


template<typename T>
void
doTask(PrintIT<Index>& P, 
       Dense<T> const& D)
    {
    auto name = std::is_same<T,Real>::value ? "Dense Real"
                                            : "Dense Cplx";
    P.printInfo(D,name,doTask(NormNoScale{},D));
     
    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.s << formatVal(P.scalefac*D.store.front()) << "\n";
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(rank);
    for(auto i : range(rank))
        gc.setRange(i,0,P.is.extent(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*D[offset(P.is,gc.i)];
        if(std::norm(val) >= Global::printScale())
            {
            P.s << "(";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                P.s << (1+gc[ii]);
                if(ii < gc.i.maxi()) P.s << ",";
                }
            P.s << ") ";

            P.s << formatVal(val) << "\n";
            }
        }
    }
template void doTask(PrintIT<Index>& P, DenseReal const& d);
template void doTask(PrintIT<Index>& P, DenseCplx const& d);

template<typename T>
Cplx
doTask(SumEls<Index>, Dense<T> const& D) 
    { 
    T sum = 0;
    for(auto& elt : D) sum += elt;
    return sum;
    }
template
Cplx
doTask(SumEls<Index>, DenseReal const& d);
template
Cplx
doTask(SumEls<Index>, DenseCplx const& d);

template<typename T1,typename T2>
void
doTask(Contract<Index> & C,
       Dense<T1> const& L,
       Dense<T2> const& R,
       ManageStore & m)
    {
    //if(not C.needresult)
    //    {
    //    m.makeNewData<ITLazy>(C.Lis,m.parg1(),C.Ris,m.parg2());
    //    return;
    //    }
    Labels Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    if(not C.Nis)
        {
        //Optimization TODO:
        //  Test different scenarios where having sortInds=true or false
        //  can improve performance. Having sorted inds can make adding
        //  quicker and let contractloop run in parallel more often in principle.
        bool sortInds = false; //whether to sort indices of result
        contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortInds);
        }
    else
        {
        Nind.resize(C.Nis.r());
        for(auto i : range(C.Nis.r()))
            {
            auto j = findindex(C.Lis,C.Nis[i]);
            if(j >= 0)
                {
                Nind[i] = Lind[j];
                }
            else
                {
                j = findindex(C.Ris,C.Nis[i]);
                Nind[i] = Rind[j];
                }
            }
        }
    auto tL = makeTenRef(L.data(),L.size(),&C.Lis);
    auto tR = makeTenRef(R.data(),R.size(),&C.Ris);
    auto rsize = area(C.Nis);
    START_TIMER(4)
    auto nd = m.makeNewData<Dense<common_type<T1,T2>>>(rsize);
    STOP_TIMER(4)
    auto tN = makeTenRef(nd->data(),nd->size(),&(C.Nis));

    START_TIMER(2)
    contract(tL,Lind,tR,Rind,tN,Nind);
    STOP_TIMER(2)

    START_TIMER(3)
    if(rsize > 1) C.scalefac = computeScalefac(*nd);
    STOP_TIMER(3)
    }
template void doTask(Contract<Index>&,DenseReal const&,DenseReal const&,ManageStore&);
template void doTask(Contract<Index>&,DenseCplx const&,DenseReal const&,ManageStore&);
template void doTask(Contract<Index>&,DenseReal const&,DenseCplx const&,ManageStore&);
template void doTask(Contract<Index>&,DenseCplx const&,DenseCplx const&,ManageStore&);

template<typename VL, typename VR>
void
doTask(NCProd<Index>& P,
       Dense<VL> const& L,
       Dense<VR> const& R,
       ManageStore& m)
    {
    Labels Lind,
          Rind,
          Nind;
    computeLabels(P.Lis,P.Lis.r(),P.Ris,P.Ris.r(),Lind,Rind);
    ncprod(P.Lis,Lind,P.Ris,Rind,P.Nis,Nind);

    auto tL = makeTenRef(L.data(),L.size(),&P.Lis);
    auto tR = makeTenRef(R.data(),R.size(),&P.Ris);
    auto rsize = area(P.Nis);
    auto nd = m.makeNewData<Dense<common_type<VL,VR>>>(rsize);
    auto tN = makeTenRef(nd->data(),nd->size(),&(P.Nis));

    ncprod(tL,Lind,tR,Rind,tN,Nind);

    if(rsize > 1) P.scalefac = computeScalefac(*nd);
    }
template void doTask(NCProd<Index>&,DenseReal const&,DenseReal const&,ManageStore&);
template void doTask(NCProd<Index>&,DenseReal const&,DenseCplx const&,ManageStore&);
template void doTask(NCProd<Index>&,DenseCplx const&,DenseReal const&,ManageStore&);
template void doTask(NCProd<Index>&,DenseCplx const&,DenseCplx const&,ManageStore&);

struct Adder
    {
    const Real f = 1.;
    Adder(Real f_) : f(f_) { }
    template<typename T1, typename T2>
    void operator()(T2 v2, T1& v1) { v1 += f*v2; }
    void operator()(Cplx v2, Real& v1) { }
    };

template<typename T1, typename T2>
void
add(PlusEQ<Index> const& P,
    Dense<T1>          & D1,
    Dense<T2>     const& D2)
    {
#ifdef DEBUG
    if(D1.size() != D2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(isTrivial(P.perm()) && std::is_same<T1,T2>::value)
        {
        auto d1 = realData(D1);
        auto d2 = realData(D2);
        daxpy_wrapper(d1.size(),P.fac(),d2.data(),1,d1.data(),1);
        }
    else
        {
        auto ref1 = makeTenRef(D1.data(),D1.size(),&P.is1());
        auto ref2 = makeTenRef(D2.data(),D2.size(),&P.is2());
        transform(permute(ref2,P.perm()),ref1,Adder{P.fac()});
        }
    }

template<typename T1, typename T2>
void
doTask(PlusEQ<Index> const& P,
       Dense<T1> const& D1,
       Dense<T2> const& D2,
       ManageStore & m)
    {
    if(isReal(D1) && isCplx(D2))
        {
        auto *ncD1 = m.makeNewData<DenseCplx>(D1.begin(),D1.end());
        add(P,*ncD1,D2);
        }
    else
        {
        auto *ncD1 = m.modifyData(D1);
        add(P,*ncD1,D2);
        }
    }
template void doTask(PlusEQ<Index> const&,Dense<Real> const&,Dense<Real> const&,ManageStore &);
template void doTask(PlusEQ<Index> const&,Dense<Real> const&,Dense<Cplx> const&,ManageStore &);
template void doTask(PlusEQ<Index> const&,Dense<Cplx> const&,Dense<Real> const&,ManageStore &);
template void doTask(PlusEQ<Index> const&,Dense<Cplx> const&,Dense<Cplx> const&,ManageStore &);

template<typename T>
void
permuteDense(Permutation const& P,
             Dense<T>    const& dA,
             IndexSet    const& Ais,
             Dense<T>         & dB,
             IndexSet    const& Bis)
    {
    auto bref = makeTenRef(dB.data(),dB.size(),&Bis);
    auto aref = makeTenRef(dA.data(),dA.size(),&Ais);
    bref &= permute(aref,P);
    }

template<typename T>
void
doTask(Order<Index> const& O,
       Dense<T> & dA)
    {
    auto dB = dA;
    permuteDense(O.perm(),dB,O.is1(),dA,O.is2());
    }
template void doTask(Order<Index> const&,Dense<Real> &);
template void doTask(Order<Index> const&,Dense<Cplx> &); 

} // namespace itensor
