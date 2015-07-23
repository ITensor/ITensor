//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/itreal.h"
#include "itensor/itdata/itdata.h"
#include "itensor/itdata/itcplx.h"
#include "itensor/itdata/itlazy.h"
#include "itensor/indexset.h"
#include "itensor/util/count.h"
#include "itensor/tensor/contract.h"
#include "itensor/matrix/lapack_wrap.h"

namespace itensor {

Cplx 
doTask(GetElt<Index> const& g, ITReal const& d)
    {
    return d[ind(g.is,g.inds)];
    }

void
doTask(SetElt<Real,Index> const& s, ITReal & d)
    {
    d[ind(s.is,s.inds)] = s.elt;
    }

void
doTask(SetElt<Cplx,Index> const& s, ITReal const& d, ManageStore & m)
    {
    auto nd = m.makeNewData<ITCplx>(d);
    nd->set(ind(s.is,s.inds),s.elt);
    }

void
doTask(FillReal const& f, ITReal & d)
    {
    std::fill(d.begin(),d.end(),f.r);
    }

void
doTask(FillCplx const& f, ITReal const& d, ManageStore & m)
    {
    m.makeNewData<ITCplx>(d.size(),f.z);
    }

void
doTask(MultCplx const& M, ITReal const& d, ManageStore & m)
    {
    auto nd = m.makeNewData<ITCplx>(d);
    (*nd) *= M.z;
    }

void
doTask(MultReal const& m, ITReal & d)
    {
    //use BLAS algorithm?
    for(auto& elt : d) elt *= m.r;
    }

Real
doTask(NormNoScale, ITReal const& d) 
    { 
    //println("In norm");
    Real nrm = 0;
    for(auto& elt : d) 
        nrm += elt*elt;
    return std::sqrt(nrm);
    }

void
doTask(Conj,ITReal const& d) { }

void
doTask(TakeReal, ITReal const& d) { }

void
doTask(TakeImag, ITReal const& d, ManageStore & m) 
    { 
    m.makeNewData<ITReal>(d.size(),0);
    }

void
doTask(PrintIT<Index>& P, 
       ITReal const& d)
    {
    P.printInfo(d,"Dense Real",doTask(NormNoScale{},d));
     
    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.printVal(P.scalefac*d.store.front());
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(0,rank-1,0);
    for(auto i : count(rank))
        gc.setInd(i,0,P.is.extent(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*d[ind(P.is,gc.i)];
        if(std::norm(val) > Global::printScale())
            {
            P.s << "(";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                P.s << (1+gc.i(ii));
                if(ii < gc.i.maxi()) P.s << ",";
                }
            P.s << ") ";

            P.printVal(val);
            }
        }
    }

Cplx
doTask(SumEls<Index>, ITReal const& d) 
    { 
    Real sum = 0;
    for(const auto& elt : d)
        sum += elt;
    return sum;
    }

void
doTask(Write & W, ITReal const& d) 
    { 
    W.writeType(StorageType::ITReal,d); 
    }

void
doTask(Contract<Index> & C,
       ITReal const& a1,
       ITReal const& a2,
       ManageStore & m)
    {
    //if(not C.needresult)
    //    {
    //    m.makeNewData<ITLazy>(C.Lis,m.parg1(),C.Ris,m.parg2());
    //    return;
    //    }
    Label Lind,
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
        for(auto i : count(C.Nis.r()))
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
    auto t1 = makeTensorRef(a1.data(),C.Lis),
         t2 = makeTensorRef(a2.data(),C.Ris);
    auto rsize = area(C.Nis);
    auto nd = m.makeNewData<ITReal>(rsize,0.);
    auto tr = makeTensorRef(nd->data(),C.Nis);
    contractloop(t1,Lind,t2,Rind,tr,Nind);

    if(rsize > 1) C.computeScalefac(*nd);
    }

void
doTask(PlusEQ<Index> const& P,
       ITReal & a1,
       ITReal const& a2)
    {
#ifdef DEBUG
    if(a1.size() != a2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(!P.hasPerm())
        {
        daxpy_wrapper(a1.size(),P.fac,a2.data(),1,a1.data(),1);
        }
    else
        {
        auto ref1 = makeTensorRef(a1.data(),P.is1());
        auto ref2 = makeTensorRef(a2.data(),P.is2());
        auto add = [f=P.fac](Real& r1, Real r2) { r1 += f*r2; };
        permute(ref2,P.perm(),ref1,add);
        }
    }

} // namespace itensor
