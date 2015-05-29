//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/itdata.h"
#include "itensor/itdata/itcplx.h"
#include "itensor/indexset.h"
#include "itensor/util/count.h"
#include "itensor/tensor/contract.h"
#include "itensor/matrix/lapack_wrap.h"

namespace itensor {

Cplx 
doTask(const GetElt<Index>& g, const ITReal& d)
    {
    return d[ind(g.is,g.inds)];
    }

void
doTask(const SetElt<Real,Index>& s, ITReal& d)
    {
    d[ind(s.is,s.inds)] = s.elt;
    }

void
doTask(const SetElt<Cplx,Index>& s, const ITReal& d, ManagePtr& mp)
    {
    auto nd = mp.makeNewData<ITCplx>(d);
    nd->set(ind(s.is,s.inds),s.elt);
    }

void
doTask(const FillReal& f, ITReal& d)
    {
    std::fill(d.begin(),d.end(),f.r);
    }

void
doTask(const FillCplx& f, const ITReal& d, ManagePtr& mp)
    {
    mp.makeNewData<ITCplx>(d.size(),f.z);
    }

void
doTask(const MultCplx& M, const ITReal& d, ManagePtr& mp)
    {
    auto nd = mp.makeNewData<ITCplx>(d);
    (*nd) *= M.z;
    }

void
doTask(const MultReal& m, ITReal& d)
    {
    //use BLAS algorithm?
    for(auto& elt : d) elt *= m.r;
    }

Real
doTask(const NormNoScale<Index>& N, const ITReal& d) 
    { 
    Real nrm = 0;
    for(auto& elt : d) 
        nrm += elt*elt;
    return std::sqrt(nrm);
    }

void
doTask(Conj,const ITReal& d) { }

void
doTask(TakeReal, const ITReal& ) { }

void
doTask(TakeImag, const ITReal& d, ManagePtr& mp) 
    { 
    mp.makeNewData<ITReal>(d.size(),0);
    }

void
doTask(PrintIT<Index>& P, const ITReal& d)
    {
    P.printInfo(d,"Dense Real",doTask(NormNoScale<Index>(P.is),d));
     
    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.printVal(P.scalefac*d.store.front());
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,P.is.dim(i)-1);

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
doTask(SumEls<Index>, const ITReal& d) 
    { 
    Real sum = 0;
    for(const auto& elt : d)
        sum += elt;
    return sum;
    }

void
doTask(Write& W, const ITReal& d) 
    { 
    W.writeType(StorageType::ITReal,d); 
    }

void
doTask(Contract<Index>& C,
       const ITReal& a1,
       const ITReal& a2,
       ManagePtr& mp)
    {
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    bool sortInds = false; //whether to sort indices of result
    contractIS(C.Lis,C.Lind,C.Ris,C.Rind,C.Nis,sortInds);
    
    Label Nind(C.Nis.r(),0);
    for(auto i : count(C.Nis.r()))
        {
        auto j = findindex(C.Lis,C.Nis[i]);
        if(j >= 0)
            {
            Nind[i] = C.Lind[j];
            }
        else
            {
            j = findindex(C.Ris,C.Nis[i]);
            Nind[i] = C.Rind[j];
            }
        }
    auto t1 = makeTensorRef(a1.data(),C.Lis),
         t2 = makeTensorRef(a2.data(),C.Ris);
    auto rsize = area(C.Nis);
    auto nd = mp.makeNewData<ITReal>(rsize,0.);
    auto tr = makeTensorRef(nd->data(),C.Nis);
    contractloop(t1,C.Lind,t2,C.Rind,tr,Nind);

    if(rsize > 1) C.computeScalefac(*nd);
    }

void
doTask(const PlusEQ<Index>& P,
       ITReal& a1,
       const ITReal& a2)
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
