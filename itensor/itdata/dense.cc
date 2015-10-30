//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/dense.h"
#include "itensor/itdata/itdata.h"
//#include "itensor/itdata/itcplx.h"
#include "itensor/itdata/itlazy.h"
#include "itensor/indexset.h"
#include "itensor/util/count.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/lapack_wrap.h"

namespace itensor {


template<typename T>
Cplx 
doTask(GetElt<Index> const& g, Dense<T> const& d)
    {
    return d[offset(g.is,g.inds)];
    }
template
Cplx 
doTask(GetElt<Index> const& g, Dense<Real> const& d);
template
Cplx 
doTask(GetElt<Index> const& g, Dense<Cplx> const& d);

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
    for(auto n : count(D.size()))
        {
        nD[n] = D[n].real();
        }
    }

void
doTask(TakeImag, DenseCplx const& D, ManageStore & m) 
    { 
    auto& nD = *m.makeNewData<DenseReal>(D.size());
    for(auto n : count(D.size()))
        {
        nD[n] = D[n].imag();
        }
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
        P.printVal(P.scalefac*D.store.front());
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(rank);
    for(auto i : count(rank))
        gc.setRange(i,0,P.is.extent(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*D[offset(P.is,gc.i)];
        if(std::norm(val) > Global::printScale())
            {
            P.s << "(";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                P.s << (1+gc[ii]);
                if(ii < gc.i.maxi()) P.s << ",";
                }
            P.s << ") ";

            P.printVal(val);
            }
        }
    }
template
void
doTask(PrintIT<Index>& P, DenseReal const& d);
template
void
doTask(PrintIT<Index>& P, DenseCplx const& d);

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

void
doTask(Write& W, DenseReal const& d)
    { 
    W.writeType(StorageType::DenseReal,d); 
    }

void
doTask(Write& W, DenseCplx const& d)
    { 
    W.writeType(StorageType::DenseCplx,d); 
    }

//template<typename T1,
//         typename T2>
//void
//doTask(Contract<Index> & C,
//       Dense<T1> const& a1,
//       Dense<T2> const& a2,
//       ManageStore & m)
//    {
//    //if(not C.needresult)
//    //    {
//    //    m.makeNewData<ITLazy>(C.Lis,m.parg1(),C.Ris,m.parg2());
//    //    return;
//    //    }
//    Label Lind,
//          Rind,
//          Nind;
//    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
//    if(not C.Nis)
//        {
//        //Optimization TODO:
//        //  Test different scenarios where having sortInds=true or false
//        //  can improve performance. Having sorted inds can make adding
//        //  quicker and let contractloop run in parallel more often in principle.
//        bool sortInds = false; //whether to sort indices of result
//        contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortInds);
//        }
//    else
//        {
//        Nind.resize(C.Nis.r());
//        for(auto i : count(C.Nis.r()))
//            {
//            auto j = findindex(C.Lis,C.Nis[i]);
//            if(j >= 0)
//                {
//                Nind[i] = Lind[j];
//                }
//            else
//                {
//                j = findindex(C.Ris,C.Nis[i]);
//                Nind[i] = Rind[j];
//                }
//            }
//        }
//    auto t1 = makeTenRef(a1.data(),a1.size(),&C.Lis),
//         t2 = makeTenRef(a2.data(),a2.size(),&C.Ris);
//    auto rsize = area(C.Nis);
//    using Ctype = common_type<decltype(t1),decltype(t2)>;
//    START_TIMER(4)
//    auto nd = m.makeNewData<Dense<Ctype>>(rsize);
//    STOP_TIMER(4)
//    auto tr = makeTenRef(nd->data(),nd->size(),&(C.Nis));
//
//    START_TIMER(2)
//    //contractloop(t1,Lind,t2,Rind,tr,Nind);
//    contract(t1,Lind,t2,Rind,tr,Nind);
//    STOP_TIMER(2)
//
//    START_TIMER(3)
//    if(rsize > 1) C.scalefac = computeScalefac(*nd);
//    STOP_TIMER(3)
//    }
//
//void
//doTask(NCProd<Index>& P,
//       Dense const& d1,
//       Dense const& d2,
//       ManageStore& m)
//    {
//    Label Lind,
//          Rind,
//          Nind;
//    computeLabels(P.Lis,P.Lis.r(),P.Ris,P.Ris.r(),Lind,Rind);
//    ncprod(P.Lis,Lind,P.Ris,Rind,P.Nis,Nind);
//
//    auto t1 = makeTenRef(d1.data(),d1.size(),&P.Lis),
//         t2 = makeTenRef(d2.data(),d2.size(),&P.Ris);
//    auto rsize = area(P.Nis);
//    auto nd = m.makeNewData<Dense>(rsize);
//    auto tr = makeTenRef(nd->data(),nd->size(),&(P.Nis));
//
//    ncprod(t1,Lind,t2,Rind,tr,Nind);
//
//    if(rsize > 1) P.scalefac = computeScalefac(*nd);
//    }

struct Adder
    {
    const Real f = 1.;
    Adder(Real f_) : f(f_) { }
    void operator()(Real v2, Real& v1) { v1 += f*v2; }
    void operator()(Real v2, Cplx& v1) { v1 += f*v2; }
    void operator()(Cplx v2, Cplx& v1) { v1 += f*v2; }
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
    if(std::is_same<T1,Real>::value && std::is_same<T2,Cplx>::value)
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
template
void
doTask(PlusEQ<Index> const& P,
       Dense<Real> const& D1,
       Dense<Real> const& D2,
       ManageStore & m);
template
void
doTask(PlusEQ<Index> const& P,
       Dense<Real> const& D1,
       Dense<Cplx> const& D2,
       ManageStore & m);
template
void
doTask(PlusEQ<Index> const& P,
       Dense<Cplx> const& D1,
       Dense<Real> const& D2,
       ManageStore & m);
template
void
doTask(PlusEQ<Index> const& P,
       Dense<Cplx> const& D1,
       Dense<Cplx> const& D2,
       ManageStore & m);

    

} // namespace itensor
