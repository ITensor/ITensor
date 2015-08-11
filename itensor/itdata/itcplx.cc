//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/itcplx.h"
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
    return d.get(offset(g.is,g.inds));
    }

void
doTask(const SetElt<Real,Index>& s, ITCplx& d)
    {
    d.set(offset(s.is,s.inds),s.elt);
    }

void
doTask(const SetElt<Cplx,Index>& s, ITCplx& d)
    {
    d.set(offset(s.is,s.inds),s.elt);
    }

void
doTask(Contract<Index>& C,
       const ITCplx& a1,
       const ITCplx& a2,
       ManageStore& m)
    {
    Label Lind,
          Rind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    bool sortInds = false; //whether to sort indices of result
    Label Nind;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortInds);

    auto rsize = area(C.Nis);
    auto nd = m.makeNewData<ITCplx>(rsize,0.);

    auto t1r = makeTenRef(a1.rstart(),C.Lis),
         t1i = makeTenRef(a1.istart(),C.Lis),
         t2r = makeTenRef(a2.rstart(),C.Ris),
         t2i = makeTenRef(a2.istart(),C.Ris);
    auto trr = makeTenRef(nd->rstart(),C.Nis),
         tri = makeTenRef(nd->istart(),C.Nis);

    contractloop(t1i,Lind,t2i,Rind,trr,Nind);
    for(auto p = nd->rstart(); p < nd->istart(); ++p) *p *= -1;
    contractloop(t1r,Lind,t2r,Rind,trr,Nind);

    contractloop(t1i,Lind,t2r,Rind,tri,Nind);
    contractloop(t1r,Lind,t2i,Rind,tri,Nind);

    if(rsize > 1) C.computeScalefac(*nd);
    }

void
realCplx(const ITReal& R,
         const IndexSet& ris,
         const ITCplx& C,
         const IndexSet& cis,
         IndexSet& Nis,
         ManageStore& m,
         Contract<Index>& Con,
         bool RealOnLeft)
    {
    Label rind,
          cind;
    computeLabels(ris,ris.r(),cis,cis.r(),rind,cind);
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    bool sortInds = false; //whether to sort indices of result
    Label Nind;
    contractIS(ris,rind,cis,cind,Nis,Nind,sortInds);

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
    

    auto nd = m.makeNewData<ITCplx>(rsize,0.);

    auto t1 = makeTenRef(R.data(),ris),
         t2r = makeTenRef(C.rstart(),cis),
         t2i = makeTenRef(C.istart(),cis);
    auto trr = makeTenRef(nd->rstart(),Nis),
         tri = makeTenRef(nd->istart(),Nis);

    contractloop(t1,rind,t2r,cind,trr,Nind);
    contractloop(t1,rind,t2i,cind,tri,Nind);

    if(rsize > 1) Con.computeScalefac(*nd);
    }


void
doTask(Contract<Index>& C,
       const ITReal& a1,
       const ITCplx& a2,
       ManageStore& m)
    {
    realCplx(a1,C.Lis,a2,C.Ris,C.Nis,m,C,true);
    }

void
doTask(Contract<Index>& C,
       const ITCplx& a1,
       const ITReal& a2,
       ManageStore& m)
    {
    realCplx(a2,C.Ris,a1,C.Lis,C.Nis,m,C,false);
    }

void
doTask(const FillReal& f, const ITCplx& d, ManageStore& m)
    {
    m.makeNewData<ITReal>(d.csize(),f.r);
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
doTask(NormNoScale, const ITCplx& d)
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
doTask(TakeReal,ITCplx& d, ManageStore& m)
    { 
    auto csize = d.csize(); //csize is half the size of d.store
    //Take resources from ITCplx to make new ITReal
    auto* nd = m.makeNewData<ITReal>(std::move(d.store));
    nd->store.resize(csize);
    }

void
doTask(TakeImag,const ITCplx& d, ManageStore& m) 
    { 
    m.makeNewData<ITReal>(d.istart(),d.iend());
    }

void
doTask(PrintIT<Index>& P, const ITCplx& d)
    {
    P.printInfo(d,"Dense Cplx",doTask(NormNoScale{},d));

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.printVal(P.scalefac*d.get(0));
        return;
        }

    if(!P.print_data) return;

    auto gc = detail::GCounter(rank);
    for(auto i : count(rank))
        gc.setRange(i,0,P.is.extent(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = P.scalefac*d.get(offset(P.is,gc.i));
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
