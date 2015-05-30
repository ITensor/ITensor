//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/itcplx.h"
//#include "itensor/itdata/itdata.h"
//#include "itensor/itdata/itreal.h"
#include "itensor/tensor/contract.h"
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/util/count.h"
#include "itensor/util/safe_ptr.h"

namespace itensor {

ITCplx::
ITCplx(const ITReal& d)
  : store(2*d.size(),0) 
    { 
    std::copy(d.begin(),d.end(),store.begin());
    }

ITCplx& ITCplx::
operator*=(const Cplx& z)
    {
    auto r = MAKE_SAFE_PTR(rstart(),csize());
    auto re = istart(); //imag start is real end
    auto i = MAKE_SAFE_PTR(istart(),csize());
    auto a = z.real(),
         b = z.imag();
    for(; r != re; ++r, ++i)
        {
        auto nr = *r*a-*i*b;
        auto ni = *i*a+*r*b;
        *r = nr;
        *i = ni;
        }
    return *this;
    }

void ITCplx::
fill(const Complex& z)
    {
    std::fill(store.begin(),store.begin()+csize(),z.real());
    std::fill(store.begin()+csize(),store.end(),z.imag());
    }

Cplx
doTask(const GetElt<Index>& g, const ITCplx& d)
    {
    return d.get(ind(g.is,g.inds));
    }

void
doTask(const SetElt<Real,Index>& s, ITCplx& d)
    {
    d.set(ind(s.is,s.inds),s.elt);
    }

void
doTask(const SetElt<Cplx,Index>& s, ITCplx& d)
    {
    d.set(ind(s.is,s.inds),s.elt);
    }

void
doTask(Contract<Index>& C,
       const ITCplx& a1,
       const ITCplx& a2,
       ManagePtr& mp)
    {
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    const bool sortInds = false; //whether to sort indices of result
    contractIS(C.Lis,C.Lind,C.Ris,C.Rind,C.Nis,sortInds);

    //Scalar tensor cases
    //if(area(C.Lis)==1)
    //    {
    //    auto nd = makeNewData<ITCplx>(a2);
    //    (*nd) *= a1.get(0);
    //    return;
    //    }
    //if(area(C.Ris)==1)
    //    {
    //    auto a1ref = modifyData(a1);
    //    a1ref *= a2.get(0);
    //    return;
    //    }
    
    Label Nind(C.Nis.r(),0);
    for(auto i : index(C.Nis))
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

    auto rsize = area(C.Nis);
    auto nd = mp.makeNewData<ITCplx>(rsize,0.);

    auto t1r = makeTensorRef(a1.rstart(),C.Lis),
         t1i = makeTensorRef(a1.istart(),C.Lis),
         t2r = makeTensorRef(a2.rstart(),C.Ris),
         t2i = makeTensorRef(a2.istart(),C.Ris);
    auto trr = makeTensorRef(nd->rstart(),C.Nis),
         tri = makeTensorRef(nd->istart(),C.Nis);

    contractloop(t1i,C.Lind,t2i,C.Rind,trr,Nind);
    for(auto p = nd->rstart(); p < nd->istart(); ++p) *p *= -1;
    contractloop(t1r,C.Lind,t2r,C.Rind,trr,Nind);

    contractloop(t1i,C.Lind,t2r,C.Rind,tri,Nind);
    contractloop(t1r,C.Lind,t2i,C.Rind,tri,Nind);

    if(rsize > 1) C.computeScalefac(*nd);
    }

void
realCplx(const ITReal& R,
         const IndexSet& ris,
         const Label& rind,
         const ITCplx& C,
         const IndexSet& cis,
         const Label& cind,
         IndexSet& Nis,
         ManagePtr& mp,
         Contract<Index>& Con,
         bool RealOnLeft)
    {
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    const bool sortInds = false; //whether to sort indices of result
    contractIS(ris,rind,cis,cind,Nis,sortInds);

    auto rsize = area(Nis);

    //if(area(ris)==1) //Real tensor is a scalar
    //    {
    //    scalefac_ = R[0];
    //    if(RealOnLeft) assignPointerRtoL();
    //    return;
    //    }
    //if(area(cis)==1) //Cplx tensor is a scalar
    //    {
    //    Error("Need to check in this case whether RealOnLeft or not");
    //    auto z = C.get(0);
    //    auto zr = z.real(),
    //         zi = z.imag();
    //    if(zi == 0)
    //        {
    //        scalefac_ = zr;
    //        assignPointerRtoL();
    //        }
    //    else
    //        {
    //        assert(size_t(rsize)==size_t(R.size()));
    //        auto nd = makeNewData<ITCplx>(rsize);
    //        auto pr = MAKE_SAFE_PTR(nd->rstart(),rsize);
    //        auto pi = MAKE_SAFE_PTR(nd->istart(),rsize);
    //        for(auto el : R)
    //            {
    //            *pr = el*zr;
    //            *pi = el*zi;
    //            ++pr;
    //            ++pi;
    //            }
    //        }
    //    return;
    //    }
    
    Label Nind(Nis.r(),0);
    for(auto i : count(Nis.r()))
        {
        auto j = findindex(ris,Nis[i]);
        if(j >= 0)
            {
            Nind[i] = rind[j];
            }
        else
            {
            j = findindex(cis,Nis[i]);
            Nind[i] = cind[j];
            }
        }

    auto nd = mp.makeNewData<ITCplx>(rsize,0.);

    auto t1 = makeTensorRef(R.data(),ris),
         t2r = makeTensorRef(C.rstart(),cis),
         t2i = makeTensorRef(C.istart(),cis);
    auto trr = makeTensorRef(nd->rstart(),Nis),
         tri = makeTensorRef(nd->istart(),Nis);

    contractloop(t1,rind,t2r,cind,trr,Nind);
    contractloop(t1,rind,t2i,cind,tri,Nind);

    if(rsize > 1) Con.computeScalefac(*nd);
    }


void
doTask(Contract<Index>& C,
       const ITReal& a1,
       const ITCplx& a2,
       ManagePtr& mp)
    {
    realCplx(a1,C.Lis,C.Lind,a2,C.Ris,C.Rind,C.Nis,mp,C,true);
    }

void
doTask(Contract<Index>& C,
       const ITCplx& a1,
       const ITReal& a2,
       ManagePtr& mp)
    {
    realCplx(a2,C.Ris,C.Rind,a1,C.Lis,C.Lind,C.Nis,mp,C,false);
    }

void
doTask(const FillReal& f, const ITCplx& d, ManagePtr& mp)
    {
    mp.makeNewData<ITReal>(d.csize(),f.r);
    }

void
doTask(const FillCplx& f, ITCplx& d)
    {
    d.fill(f.z);
    }

void
doTask(const MultCplx& M, ITCplx& d) 
    { 
    d *= M.z; 
    }

void
doTask(const MultReal& m, ITCplx& d)
    {
    for(auto& elt : d) elt *= m.r;
    }

Real
doTask(const NormNoScale<Index>& N, const ITCplx& d)
    { 
    Real nrm = 0;
    for(auto& elt : d)
        nrm += std::norm(elt); //conj(elt)*elt
    return std::sqrt(nrm);
    }

void
doTask(Conj, ITCplx& d) 
    { 
    auto i = MAKE_SAFE_PTR(d.istart(),d.csize());
    auto* ie = d.iend();
    for(; i != ie; ++i) *i *= -1;
    }

void
doTask(TakeReal,ITCplx& d, ManagePtr& mp)
    { 
    auto csize = d.csize(); //csize is half the size of d.store
    //Take resources from ITCplx to make new ITReal
    auto* nd = mp.makeNewData<ITReal>(std::move(d.store));
    nd->store.resize(csize);
    }

void
doTask(TakeImag,const ITCplx& d, ManagePtr& mp) 
    { 
    mp.makeNewData<ITReal>(d.istart(),d.iend());
    }

void
doTask(PrintIT<Index>& P, const ITCplx& d)
    {
    P.printInfo(d,"Dense Cplx",doTask(NormNoScale<Index>(P.is),d));

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.printVal(P.scalefac*d.get(0));
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,P.is.dim(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*d.get(ind(P.is,gc.i));
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

bool
doTask(CheckComplex,const ITCplx& d) { return true; }

Cplx
doTask(SumEls<Index>, const ITCplx& d) 
    { 
    Real rsum = 0,
         isum = 0;
    auto p = MAKE_SAFE_PTR(d.rstart(),d.csize());
    for(; p != d.istart(); ++p) rsum += *p;
    p = MAKE_SAFE_PTR(d.istart(),d.csize());
    for(; p != d.iend(); ++p)   isum += *p;
    return Cplx(rsum,isum);
    }

void
doTask(Write& W, const ITCplx& d) 
    { 
    W.writeType(StorageType::ITCplx,d); 
    }

} // namespace itensor
