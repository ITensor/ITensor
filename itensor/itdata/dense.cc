//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#include "itensor/itdata/dense.h"
#include "itensor/itdata/itdata.h"
//#include "itensor/itdata/itlazy.h"
#include "itensor/indexset.h"
#include "itensor/util/iterate.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/util/tensorstats.h"

using std::move;
using std::string;
using std::vector;

namespace itensor {

const char*
typeNameOf(DenseReal const& d) { return "DenseReal"; }
const char*
typeNameOf(DenseCplx const& d) { return "DenseCplx"; }

Cplx 
doTask(GetElt const& g, DenseReal const& d)
    {
    return Cplx(d[offset(g.is,g.inds)],0.);
    }

Cplx 
doTask(GetElt const& g, DenseCplx const& d)
    {
    return d[offset(g.is,g.inds)];
    }

template<typename E, typename T>
struct SetEltHelper
    {
    void static
    set(SetElt<E> const& S, Dense<T> const& D, ManageStore& m)
        {
        auto& Dnc = *m.modifyData(D);
        Dnc[offset(S.is,S.inds)] = S.elt;
        }
    };
template<>
struct SetEltHelper<Cplx,Real>
    {
    void static
    set(SetElt<Cplx> const& S, DenseReal const& D, ManageStore & m)
        {
        auto& nd = *m.makeNewData<DenseCplx>(D.begin(),D.end());
        nd[offset(S.is,S.inds)] = S.elt;
        }
    };

template<typename E, typename T>
void
doTask(SetElt<E> const& S, Dense<T> const& D, ManageStore & m)
    {
    SetEltHelper<E,T>::set(S,D,m);
    }
template
void
doTask(SetElt<Real> const& S, DenseReal const& D, ManageStore & m);
template
void
doTask(SetElt<Real> const& S, DenseCplx const& D, ManageStore & m);
template
void
doTask(SetElt<Cplx> const& S, DenseReal const& D, ManageStore & m);
template
void
doTask(SetElt<Cplx> const& S, DenseCplx const& D, ManageStore & m);


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

void
doTask(MakeCplx const&, Dense<Cplx> & D)
    {
    //do nothing, already complex
    }
void
doTask(MakeCplx const&, Dense<Real> const& D, ManageStore & m)
    {
    //convert data to complex
    m.makeNewData<DenseCplx>(D.begin(),D.end());
    }

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


template<typename T>
void
doTask(PrintIT& P, 
       Dense<T> const& D)
    {
    auto name = std::is_same<T,Real>::value ? "Dense Real"
                                            : "Dense Cplx";
    P.printInfo(D,name,doTask(NormNoScale{},D));
     
    auto ord = P.is.order();
    if(ord == 0) 
        {
        P.s << "  ";
        P.s << formatVal(P.scalefac*D.store.front()) << "\n";
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(ord);
    for(auto i : range(ord))
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
template void doTask(PrintIT& P, DenseReal const& d);
template void doTask(PrintIT& P, DenseCplx const& d);

template<typename T>
Cplx
doTask(SumEls, Dense<T> const& D) 
    { 
    T sum = 0;
    for(auto& elt : D) sum += elt;
    return sum;
    }
template
Cplx
doTask(SumEls, DenseReal const& d);
template
Cplx
doTask(SumEls, DenseCplx const& d);

template<typename T1,typename T2>
void
doTask(Contract & C,
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
    computeLabels(C.Lis,C.Lis.order(),C.Ris,C.Ris.order(),Lind,Rind);
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
        Nind.resize(C.Nis.order());
        for(auto i : range(C.Nis.order()))
            {
            auto j = indexPosition(C.Lis,C.Nis[i]);
            if(j >= 0)
                {
                Nind[i] = Lind[j];
                }
            else
                {
                j = indexPosition(C.Ris,C.Nis[i]);
                Nind[i] = Rind[j];
                }
            }
        }
    auto tL = makeTenRef(L.data(),L.size(),&C.Lis);
    auto tR = makeTenRef(R.data(),R.size(),&C.Ris);
    auto rsize = dim(C.Nis);
TIMER_START(40);
    // Create a Dense storage with undefined data, since it will be
    // overwritten anyway
    auto nd = m.makeNewData<Dense<common_type<T1,T2>>>(undef,rsize);
TIMER_STOP(40);
    auto tN = makeTenRef(nd->data(),nd->size(),&(C.Nis));

#ifdef COLLECT_TSTATS
    tstats(tL,Lind,tR,Rind,tN,Nind);
#endif

START_TIMER(41);
    contract(tL,Lind,tR,Rind,tN,Nind);
STOP_TIMER(41);


#ifdef USESCALE
    if(rsize > 1) C.scalefac = computeScalefac(*nd);
#endif
    }
template void doTask(Contract&,DenseReal const&,DenseReal const&,ManageStore&);
template void doTask(Contract&,DenseCplx const&,DenseReal const&,ManageStore&);
template void doTask(Contract&,DenseReal const&,DenseCplx const&,ManageStore&);
template void doTask(Contract&,DenseCplx const&,DenseCplx const&,ManageStore&);

template<typename VL, typename VR>
void
doTask(NCProd& P,
       Dense<VL> const& L,
       Dense<VR> const& R,
       ManageStore& m)
    {
    Labels Lind,
          Rind,
          Nind;
    computeLabels(P.Lis,P.Lis.order(),P.Ris,P.Ris.order(),Lind,Rind);
    ncprod(P.Lis,Lind,P.Ris,Rind,P.Nis,Nind);

    auto tL = makeTenRef(L.data(),L.size(),&P.Lis);
    auto tR = makeTenRef(R.data(),R.size(),&P.Ris);
    auto rsize = dim(P.Nis);
    auto nd = m.makeNewData<Dense<common_type<VL,VR>>>(rsize);
    auto tN = makeTenRef(nd->data(),nd->size(),&(P.Nis));

    ncprod(tL,Lind,tR,Rind,tN,Nind);

#ifdef USESCALE
    if(rsize > 1) P.scalefac = computeScalefac(*nd);
#endif
    }
template void doTask(NCProd&,DenseReal const&,DenseReal const&,ManageStore&);
template void doTask(NCProd&,DenseReal const&,DenseCplx const&,ManageStore&);
template void doTask(NCProd&,DenseCplx const&,DenseReal const&,ManageStore&);
template void doTask(NCProd&,DenseCplx const&,DenseCplx const&,ManageStore&);

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
add(PlusEQ const& P,
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
        daxpy_wrapper(d1.size(),P.alpha(),d2.data(),1,d1.data(),1);
        }
    else
        {
        auto ref1 = makeTenRef(D1.data(),D1.size(),&P.is1());
        auto ref2 = makeTenRef(D2.data(),D2.size(),&P.is2());
        transform(permute(ref2,P.perm()),ref1,Adder{P.alpha()});
        }
    }

template<typename T1, typename T2>
void
doTask(PlusEQ const& P,
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
template void doTask(PlusEQ const&,Dense<Real> const&,Dense<Real> const&,ManageStore &);
template void doTask(PlusEQ const&,Dense<Real> const&,Dense<Cplx> const&,ManageStore &);
template void doTask(PlusEQ const&,Dense<Cplx> const&,Dense<Real> const&,ManageStore &);
template void doTask(PlusEQ const&,Dense<Cplx> const&,Dense<Cplx> const&,ManageStore &);

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
doTask(Order const& O,
       Dense<T> & dA)
    {
    auto dB = dA;
    permuteDense(O.perm(),dB,O.is1(),dA,O.is2());
    }
template void doTask(Order const&,Dense<Real> &);
template void doTask(Order const&,Dense<Cplx> &); 

#ifdef ITENSOR_USE_HDF5

const char*
juliaTypeNameOf(DenseReal const& d) { return "Dense{Float64}"; }
const char*
juliaTypeNameOf(DenseCplx const& d) { return "Dense{ComplexF64}"; }

template<typename V>
void
h5_write(h5::group parent, std::string const& name, Dense<V> const& D)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type",juliaTypeNameOf(D),true);
    h5_write_attribute(g,"version",long(1));
    auto data = std::vector<V>(D.store.begin(),D.store.end());
    h5_write(g,"data",data);
    }
template void h5_write(h5::group, std::string const&, Dense<Real> const& D);
template void h5_write(h5::group, std::string const&, Dense<Cplx> const& D);

template<typename V>
void
h5_read(h5::group parent, std::string const& name, Dense<V> & D)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<string>(g,"type");
    if(type != juliaTypeNameOf(D)) Error(format("Group does not contain %s data in HDF5 file",typeNameOf(D)));
    auto data = h5_read<vector<V>>(g,"data");
    D = Dense<V>(move(data));
    }
template void h5_read(h5::group parent, std::string const& name, Dense<Real> & D);
template void h5_read(h5::group parent, std::string const& name, Dense<Cplx> & D);


#endif //ITENSOR_USE_HDF5

} // namespace itensor
