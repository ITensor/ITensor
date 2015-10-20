//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/itcplx.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"

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
       ITCplx const& a1,
       ITCplx const& a2,
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

    auto t1r = makeTenRef(a1.rstart(),a1.csize(),&C.Lis),
         t1i = makeTenRef(a1.istart(),a1.csize(),&C.Lis),
         t2r = makeTenRef(a2.rstart(),a2.csize(),&C.Ris),
         t2i = makeTenRef(a2.istart(),a2.csize(),&C.Ris);
    auto trr = makeTenRef(nd->rstart(),nd->csize(),&C.Nis),
         tri = makeTenRef(nd->istart(),nd->csize(),&C.Nis);

    contract(t1i,Lind,t2i,Rind,trr,Nind,1,1);
    for(auto p = nd->rstart(); p < nd->istart(); ++p) *p *= -1;
    contract(t1r,Lind,t2r,Rind,trr,Nind,1,1);

    contract(t1i,Lind,t2r,Rind,tri,Nind,1,1);
    contract(t1r,Lind,t2i,Rind,tri,Nind,1,1);

    //contractloop(t1i,Lind,t2i,Rind,trr,Nind);
    //for(auto p = nd->rstart(); p < nd->istart(); ++p) *p *= -1;
    //contractloop(t1r,Lind,t2r,Rind,trr,Nind);
    //contractloop(t1i,Lind,t2r,Rind,tri,Nind);
    //contractloop(t1r,Lind,t2i,Rind,tri,Nind);

    if(rsize > 1) C.computeScalefac(*nd);
    }

void
realCplx(ITReal const& R,
         IndexSet const& ris,
         ITCplx const& C,
         IndexSet const& cis,
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
    auto nd = m.makeNewData<ITCplx>(rsize,0.);

    auto t1 = makeTenRef(R.data(),R.size(),&ris),
         t2r = makeTenRef(C.rstart(),C.csize(),&cis),
         t2i = makeTenRef(C.istart(),C.csize(),&cis);
    auto trr = makeTenRef(nd->rstart(),nd->csize(),&Nis),
         tri = makeTenRef(nd->istart(),nd->csize(),&Nis);

    contract(t1,rind,t2r,cind,trr,Nind);
    contract(t1,rind,t2i,cind,tri,Nind);

    //contractloop(t1,rind,t2r,cind,trr,Nind);
    //contractloop(t1,rind,t2i,cind,tri,Nind);

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
doTask(PlusEQ<Index> const& P,
       ITCplx & a1,
       ITCplx const& a2)
    {
    if(isTrivial(P.perm()))
        {
        daxpy_wrapper(a1.size(),P.fac(),a2.data(),1,a1.data(),1);
        }
    else
        {
        auto t1r = makeTenRef(a1.rstart(),a1.csize(),&P.is1());
        auto t1i = makeTenRef(a1.istart(),a1.csize(),&P.is1());
        auto t2r = makeTenRef(a2.rstart(),a2.csize(),&P.is2());
        auto t2i = makeTenRef(a2.istart(),a2.csize(),&P.is2());
        auto f = P.fac();
        auto add = [f](Real r2, Real& r1) { r1 += f*r2; };
        transform(permute(t2r,P.perm()),t1r,add);
        transform(permute(t2i,P.perm()),t1i,add);
        }
    }

//void
//addCplxReal(ITCplx             & C,
//            IndexSet      const& cis,
//            ITReal        const& R,
//            IndexSet      const& ris,
//            PlusEQ<Index> const& P,
//            bool Cleft)
//    {
//    }

void
doTask(PlusEQ<Index> & P,
       ITReal const& a1,
       ITCplx const& a2,
       ManageStore& m)
    {
    auto *nd = m.makeNewData<ITCplx>(a1);
    doTask(P,*nd,a2);

    //auto *nd = m.makeNewData<ITCplx>(a2);
    //if(P.fac() != 1.0)
    //    {
    //    doTask(MultReal{P.fac},*nd);
    //    P.fac() = 1.0;
    //    }
    //addCplxReal(*nd,P.is2(),a1,P.is1(),P,false);
    //P.switchIndSet = true;
    }

void
doTask(PlusEQ<Index> const& P,
       ITCplx & C,
       ITReal const& R)
    {
    if(isTrivial(P.perm()))
        {
        daxpy_wrapper(C.csize(),P.fac(),R.data(),1,C.rstart(),1);
        }
    else
        {
        auto Cr = makeTenRef(C.rstart(),C.csize(),&P.is1());
        auto Rr = makeTenRef(R.data(),R.size(),&P.is2());
        auto f = P.fac();
        auto add = [f](Real r2, Real& r1) { r1 += f*r2; };
        transform(permute(Rr,P.perm()),Cr,add);
        }
    }

void
doTask(FillReal const& f, ITCplx const& d, ManageStore& m)
    {
    m.makeNewData<ITReal>(d.csize(),f.r);
    }

void
doTask(const FillCplx& f, ITCplx& d)
    {
    d.fill(f.z);
    }

void
doTask(MultCplx const& M, ITCplx& d) 
    { 
    d *= M.z; 
    }

void
doTask(MultReal const& m, ITCplx& d)
    {
    for(auto& elt : d) elt *= m.r;
    }

Real
doTask(NormNoScale, ITCplx const& d)
    { 
    return dnrm2_wrapper(d.size(),d.data());
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
