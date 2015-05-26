//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/detail/printing.h"
#include "itensor/itensor.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"
#include "itensor/util/safe_ptr.h"

using std::array;
using std::ostream;
using std::vector;
using std::make_shared;

namespace itensor {

//
// ITensor Constructors
//


template<>
ITensor::
ITensorT(const Index& i1) 
    :
    is_(i1),
    store_(std::make_shared<ITDataType<ITReal>>(i1.m(),0.)),
    scale_(1.)
	{ }


template<>
ITensor::
ITensorT(const Index& i1,const Index& i2) 
    :
    is_(i1,i2),
    store_(std::make_shared<ITDataType<ITReal>>(i1.m()*i2.m(),0.)),
    scale_(1.)
	{ }
    
template<>
ITensor::
ITensorT(Cplx val) 
    :
    scale_(1.)
    { 
    if(val.imag() == 0)
        store_ = std::make_shared<ITDataType<ITReal>>(1,val.real());
    else
        store_ = std::make_shared<ITDataType<ITCplx>>(1,val);
    //if(val.imag() == 0)
    //    store_ = std::make_shared<ITDataType<ITDiag<Real>>>(val.real());
    //else
    //    store_ = std::make_shared<ITDataType<ITDiag<Cplx>>>(val);
    }

template<>
ITensor::
ITensorT(IndexSet iset,
        storage_ptr&& pdat,
        const LogNumber& scale)
    :
    is_(std::move(iset)),
    store_(std::move(pdat)),
    scale_(scale)
    {
    }

vector<Index>
computeNewInds(const IndexSet& Lis,
               const Label& Lind,
               const IndexSet& Ris,
               const Label& Rind,
               size_t size_hint = 0)
    {
    vector<Index> newind;
    if(size_hint > 0) newind.reserve(size_hint);
    for(int j = 0; j < Lis.r(); ++j)
        {
        if(Lind[j] > 0) newind.push_back(Lis[j]);
        }
    for(int j = 0; j < Ris.r(); ++j)
        {
        if(Rind[j] > 0) newind.push_back(Ris[j]);
        }
    return newind;
    }


template<typename Data>
Real
computeScalefac(Data& dat)
    {
    Real scalefac = 0;
    for(auto elt : dat) scalefac += elt*elt;
    scalefac = std::sqrt(scalefac);
    if(scalefac == 0) return 0;
    for(auto& elt : dat) elt /= scalefac;
    return scalefac;
    }
  
void
doTask(Contract& C,
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

    if(rsize > 1) C.scalefac = computeScalefac(*nd);
    }


Real
realCplx(const ITReal& R,
         const IndexSet& ris,
         const Label& rind,
         const ITCplx& C,
         const IndexSet& cis,
         const Label& cind,
         IndexSet& Nis,
         ManagePtr& mp,
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

    if(rsize > 1) return computeScalefac(*nd);
    return NAN;
    }


void
doTask(Contract& C,
       const ITReal& a1,
       const ITCplx& a2,
       ManagePtr& mp)
    {
    C.scalefac = realCplx(a1,C.Lis,C.Lind,a2,C.Ris,C.Rind,C.Nis,mp,true);
    }

void
doTask(Contract& C,
       const ITCplx& a1,
       const ITReal& a2,
       ManagePtr& mp)
    {
    C.scalefac = realCplx(a2,C.Ris,C.Rind,a1,C.Lis,C.Lind,C.Nis,mp,false);
    }


void
doTask(Contract& C,
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

    //Scalar tensor cases
    //if(area(C.Lis)==1)
    //    {
    //    scalefac_ = a1[0];
    //    assignPointerRtoL();
    //    return;
    //    }
    //if(area(C.Ris)==1)
    //    {
    //    scalefac_ = a2[0];
    //    return;
    //    }
    
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

    //PRI(C.Lind);
    //PRI(C.Rind);
    //PRI(Nind);

    auto t1 = makeTensorRef(a1.data(),C.Lis),
         t2 = makeTensorRef(a2.data(),C.Ris);
    auto rsize = area(C.Nis);
    auto nd = mp.makeNewData<ITReal>(rsize,0.);
    auto tr = makeTensorRef(nd->data(),C.Nis);
    contractloop(t1,C.Lind,t2,C.Rind,tr,Nind);

    if(rsize > 1) C.scalefac = computeScalefac(*nd);
    }

void
diagDense(const ITDiag<Real>& d,
          const IndexSet& dis,
          const Label& dind,
          const ITReal& t,
          const IndexSet& tis,
          const Label& tind,
          IndexSet& Nis,
          ManagePtr& mp,
          bool RealonLeft)
    {
    //if(area(tis)==1) //Dense tensor is a scalar
    //    {
    //    scalefac_ = R[0];
    //    if(RealOnLeft) assignPointerRtoL();
    //    return;
    //    }
    //if(area(dis)==1) //Diag tensor is a scalar
    //    {
    //    Error("Not implemented");
    //    return;
    //    }

    long t_cstride = 0; //total t-stride of contracted inds of t
    size_t ntu = 0; //number uncontracted inds of t
    assert(tind.size() == tis.size());
    for(auto j : index(tind))
        {
        //if index j is contracted, add its stride to t_cstride:
        if(tind[j] < 0) t_cstride += tis.stride(j);
        else            ++ntu;
        }

    long d_ustride = 0; //total result-stride of uncontracted inds of d
    for(auto i : index(Nis))
        {
        auto j = findindex(dis,Nis[i]);
        if(j >= 0) d_ustride += Nis.stride(i);
        }

    auto dsize = size_t(minM(dis));

    if(ntu > 0)
        {
        vector<long> tstride(ntu,0),
                     rstride(ntu,0);
        detail::GCounter C(0,ntu,0);
        size_t n = 0;
        for(auto j : index(tind))
            {
            if(tind[j] > 0)
                {
#ifdef DEBUG
                if(n >= ntu) Error("n out of range");
#endif
                C.setInd(n,0,tis.dim(j)-1);
                tstride.at(n) = tis.stride(j);
                auto k = findindex(Nis,tis[j]);
#ifdef DEBUG
                if(k < 0) Error("Index not found");
#endif
                rstride.at(n) = Nis.stride(k);
                ++n;
                }
            }
        auto nd = mp.makeNewData<ITReal>(area(Nis),0.);
        auto pr = MAKE_SAFE_PTR(nd->data(),nd->size());
        auto pt = MAKE_SAFE_PTR(t.data(),t.size());

        if(d.allSame())
            {
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(auto i : count(ntu))
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(auto J : count(dsize))
                    {
                    pr[J*d_ustride+roffset] += d.val*pt[J*t_cstride+toffset];
                    }
                }
            }
        else
            {
            auto pd = MAKE_SAFE_PTR(d.data(),d.size());
            assert(d.size() == dsize);
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(auto i : count(ntu))
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(auto J : count(dsize))
                    {
                    pr[J*d_ustride+roffset] += pd[J]*pt[J*t_cstride+toffset];
                    }
                }
            }
        }
    else
        {
        //all of t's indices contracted with d
        //result will be diagonal
        if(d_ustride == 0) //all of d's inds contracted
            {
            // o scalar if all of d's inds contracted also
            Real val = 0;
            auto pt = MAKE_SAFE_PTR(t.data(),t.size());
            if(d.allSame())
                {
                for(auto J : count(dsize))
                    val += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(dsize == d.size());
                auto pd = MAKE_SAFE_PTR(d.data(),d.size());
                for(auto J : count(dsize))
                    val += pd[J]*pt[J*t_cstride];
                }
            mp.makeNewData<ITDiag<Real>>(val);
            }
        else //some of d's inds uncontracted
            {
            // o element-wise product of d's data and t's diagonal
            auto nd = mp.makeNewData<ITDiag<Real>>(dsize,0.);
            auto pr = MAKE_SAFE_PTR(nd->data(),nd->size());
            auto pt = MAKE_SAFE_PTR(t.data(),t.size());
            if(d.allSame())
                {
                for(auto J : count(dsize))
                    pr[J] += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(dsize == d.size());
                auto pd = MAKE_SAFE_PTR(d.data(),d.size());
                for(auto J : count(dsize))
                    pr[J] += pd[J]*pt[J*t_cstride];
                }
            }
        }
    }

void
doTask(Contract& C,
       const ITReal& t,
       const ITDiag<Real>& d,
       ManagePtr& mp)
    { 
    contractIS(C.Lis,C.Lind,C.Ris,C.Rind,C.Nis,true);
    diagDense(d,C.Ris,C.Rind,t,C.Lis,C.Lind,C.Nis,mp,true);
    }
void
doTask(Contract& C,
       const ITDiag<Real>& d,
       const ITReal& t,
       ManagePtr& mp)
    {
    contractIS(C.Lis,C.Lind,C.Ris,C.Rind,C.Nis,true);
    diagDense(d,C.Lis,C.Lind,t,C.Ris,C.Rind,C.Nis,mp,false);
    }

void
combine(const ITReal& d,
        const IndexSet& dis,
        const IndexSet& Cis,
        IndexSet& Nis,
        ManagePtr& mp)
    {
    //TODO: try to make use of Lind,Rind label vectors
    //      to simplify combine logic
    const auto& cind = Cis[0];
    auto jc = findindex(dis,cind);
    if(jc >= 0) //has cind
        {
        //dis has cind, replace with other inds
        vector<Index> newind;
        newind.reserve(dis.r()+Cis.r()-2);
        for(auto j : count(dis.r()))
            if(j == jc)
                {
                for(size_t k = 1; k < Cis.size(); ++k)
                    newind.push_back(Cis[k]);
                }
            else
                {
                newind.push_back(dis[j]);
                }
        Nis = IndexSet(move(newind));
        }
    else
        {
        //dis doesn't have cind, replace
        //Cis[1], Cis[2], ... with cind
        //may need to permute
        auto J1 = findindex(dis,Cis[1]);
        if(J1 < 0) 
            {
            println("IndexSet of dense tensor = \n",dis);
            println("IndexSet of combiner/delta = \n",Cis);
            Error("No contracted indices in combiner-tensor product");
            }
        //Check if Cis[1],Cis[2],... are grouped together (contiguous)
        //and in same order as on combiner
        bool contig_sameord = true;
        int c = 2;
        for(int j = J1+1; c < Cis.r() && j < dis.r(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                contig_sameord = false;
                break;
                }
        if(c != Cis.r()) contig_sameord = false;

        //printfln("%s:",contig_sameord?"Contig":"Not Contig");
        //println("  dis = ",dis);
        //println("  Cis = ",Cis);

        if(contig_sameord)
            {
            vector<Index> newind;
            newind.reserve(dis.r()-Cis.r()+2);
            for(int j = 0; j < J1; ++j) 
                newind.push_back(dis[j]);
            newind.push_back(cind);
            for(int j = J1+Cis.r()-1; j < dis.r(); ++j) 
                newind.push_back(dis[j]);
            assert(newind.size() == size_t(dis.r()-Cis.r()+2));
            Nis = IndexSet(move(newind));
            }
        else
            {
            Permutation P(dis.r());
            //Set P destination values to -1 to mark
            //indices that need to be assigned destinations:
            for(int i = 0; i < P.size(); ++i) P.setFromTo(i,-1);

            //permute combined indices to the front, in same
            //order as in Cis:
            long ni = 0;
            for(auto c : count(1,Cis.r()))
                {
                auto j = findindex(dis,Cis[c]);
                if(j < 0) 
                    {
                    println("IndexSet of dense tensor =\n  ",dis);
                    println("IndexSet of combiner/delta =\n  ",Cis);
                    println("Missing index: ",Cis[c]);
                    Error("Combiner: missing index");
                    }
                P.setFromTo(j,ni++);
                }
            //permute uncombined indices to back, keeping relative order:
            Range::storage_type pdims(dis.r());
            vector<Index> newind;
            newind.reserve(dis.r()-Cis.r()+2);
            newind.push_back(cind);
            for(auto j : count(dis.r()))
                {
                if(P.dest(j) == -1) 
                    {
                    P.setFromTo(j,ni++);
                    newind.push_back(dis[j]);
                    }
#ifdef DEBUG
                pdims.at(P.dest(j)).dim = dis[j].m();
#else
                pdims[P.dest(j)].dim = dis[j].m();
#endif
                }
            assert(newind.size()==size_t(dis.r()-Cis.r()+2));
            Range rr(move(pdims));
            Nis = IndexSet(move(newind));
            auto nd = mp.makeNewData<ITReal>(area(Nis));
            auto tr = makeTensorRef(nd->data(),rr);
            auto td = makeTensorRef(d.data(),dis);
            permute(td,P,tr);
            }
        }
    }

void
doTask(Contract& C,
       const ITReal& d,
       const ITCombiner& cmb,
       ManagePtr& mp)
    {
    combine(d,C.Lis,C.Ris,C.Nis,mp);
    }
void
doTask(Contract& C,
       const ITCombiner& cmb,
       const ITReal& d,
       ManagePtr& mp)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,mp);
    if(!mp.newData()) mp.assignPointerRtoL();
    }


ITensor&
operator*=(ITensor& A, const ITensor& B)
    {
    if(!A || !B)
        Error("Default constructed ITensor in product");

    if(&A == &B)
        {
        A = ITensor(sqr(norm(A)));
        return A;
        }

    auto& Lis = A.inds();
    auto& Ris = B.inds();

    Label Lind,
          Rind;
    //auto ncont =
    computeLabels(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind);

    //Check if other is a scalar (modulo m==1 inds)
    //if(Ris.rn() == 0)
    //    {
    //    auto nuniq = Lis.r()+Ris.r()-2*ncont;
    //    operator*=(other.cplx());
    //    is_ = IndexSet(computeNewInds(Lis,Lind,Ris,Rind,nuniq));
    //    return *this;
    //    }
    //Check if this is a scalar (modulo m==1 inds)
    //---> This case is problematic to implement this way.
    //     For example, what if other has ITCombiner storage?
    //if(Lis.rn() == 0)
    //    {
    //    auto nuniq = Lis.r()+Ris.r()-2*ncont;
    //    auto newind = computeNewInds(Lis,Lind,Ris,Rind,nuniq);
    //    operator=(other*cplx());
    //    is_ = IndexSet(std::move(newind));
    //    return *this;
    //    }

    auto nstore = A.store();
    auto C = doTask(Contract{Lis,Lind,Ris,Rind},nstore,B.store());

    auto nscale = A.scale() * B.scale();
    if(!std::isnan(C.scalefac)) nscale *= C.scalefac;

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    A = ITensor(C.Nis,std::move(nstore),nscale);

    return A;
    }

void
doTask(const FillReal& f, ITReal& d)
    {
    std::fill(d.begin(),d.end(),f.r);
    }

void
doTask(const FillReal& f, const ITCplx& d, ManagePtr& mp)
    {
    mp.makeNewData<ITReal>(d.csize(),f.r);
    }

template<typename T>
void
doTask(const FillReal& f, const ITDiag<T>& d, ManagePtr& mp)
    {
    mp.makeNewData<ITDiag<Real>>(f.r);
    }

void
doTask(const FillCplx& f, const ITReal& d, ManagePtr& mp)
    {
    mp.makeNewData<ITCplx>(d.size(),f.z);
    }

void
doTask(const FillCplx& f, ITCplx& d)
    {
    d.fill(f.z);
    }

template<typename T>
void
doTask(const FillCplx& f, const ITDiag<T>& d, ManagePtr& mp)
    {
    mp.makeNewData<ITDiag<Cplx>>(f.z);
    }

template<>
ITensor& ITensor::
fill(Cplx z)
    {
    if(!bool(*this)) return *this;
    scale_ = LogNumber(1.);
    if(z.imag() == 0)
        doTask(FillReal{z.real()},store_);
    else
        doTask(FillCplx{z},store_);
    return *this;
    }


void
doTask(const MultCplx& M, const ITReal& d, ManagePtr& mp)
    {
    auto nd = mp.makeNewData<ITCplx>(d);
    (*nd) *= M.z;
    }

void
doTask(const MultCplx& M, ITCplx& d) 
    { 
    d *= M.z; 
    }

ITensor& 
operator*=(ITensor& T, Cplx z)
    {
    if(z.imag() == 0) return operator*=(T,z.real());
    doTask(MultCplx{z},T.store());
    return T;
    }


template<>
void ITensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    doTask(MultReal{scale_.real0()},store_);
    scale_ = newscale;
    }

void
plusEqData(Real fac, Real *d1, const Real *d2, LAPACK_INT size)
    {
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,d2,&inc,d1,&inc);
    }


void
doTask(const PlusEQ& P,
       ITReal& a1,
       const ITReal& a2)
    {
#ifdef DEBUG
    if(a1.size() != a2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(!P.hasPerm())
        {
        plusEqData(P.fac,a1.data(),a2.data(),a1.size());
        }
    else
        {
        auto ref1 = makeTensorRef(a1.data(),P.is1());
        auto ref2 = makeTensorRef(a2.data(),P.is2());
        auto add = [f=P.fac](Real& r1, Real r2) { r1 += f*r2; };
        permute(ref2,P.perm(),ref1,add);
        }
    }

void
doTask(const PlusEQ& P,
       ITDiag<Real>& a1,
       const ITDiag<Real>& a2)
    {
#ifdef DEBUG
    if(a1.size() != a2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(a1.allSame() || a2.allSame()) Error("ITDiag plusEq allSame case not implemented");
    plusEqData(P.fac,a1.data(),a2.data(),a1.size());
    }

ITensor&
operator+=(ITensor& A, const ITensor& B)
    {
    if(!A) Error("Calling += on default constructed ITensor");
    if(!B) Error("Right-hand-side of ITensor += is default constructed");
    if(&A == &B) return operator*=(A,2.);
    if(A.scale().isZero()) return A.operator=(B);

    PlusEQ::permutation P(A.inds().size());
#ifdef DEBUG
    try {
        calc_permutation(B.inds(),A.inds(),P);
        }
    catch(const std::exception& e)
        {
        Print(A);
        Print(B);
        Error("ITensor::operator+=: different Index structure");
        }
#else
    calc_permutation(B.inds(),A.inds(),P);
#endif

    Real scalefac = 1;
    if(A.scale().magnitudeLessThan(B.scale())) 
        {
        A.scaleTo(B.scale()); 
        }
    else
        {
        scalefac = (B.scale()/A.scale()).real();
        }

    if(isTrivial(P))
        {
        doTask(PlusEQ{scalefac},A.store(),B.store());
        }
    else
        {
        doTask(PlusEQ{P,A.inds(),B.inds(),scalefac},A.store(),B.store());
        }

    return A;
    } 

ITensor&
operator-=(ITensor& A, const ITensor& B)
    {
    if(&A == &B) 
        { 
        A.scale() = 0; 
        A.fill(0);
        return A;
        }
    A.scale().negate();
    operator+=(A,B); 
    A.scale().negate();
    return A;
    }



void
doTask(const MultReal& m, ITReal& d)
    {
    //use BLAS algorithm?
    for(auto& elt : d) elt *= m.r;
    }

void
doTask(const MultReal& m, ITCplx& d)
    {
    for(auto& elt : d) elt *= m.r;
    }

template<typename T>
void
doTask(const MultReal& m, ITDiag<T>& d)
    {
    d.val *= m.r;
    //use BLAS algorithm?
    for(auto& elt : d.store) elt *= m.r;
    }



template<typename Container>
Real
vec_norm(const Container& v)
    {
    Real nrm = 0;
    for(const auto& elt : v)
        nrm += std::norm(elt); //conj(elt)*elt
    return std::sqrt(nrm);
    }

Real
doTask(const NormNoScale& N, const ITReal& d) { return vec_norm(d); }

Real
doTask(const NormNoScale& N, const ITCplx& d) { return vec_norm(d); }

template<typename T>
Real
doTask(const NormNoScale& N, const ITDiag<T>& d)
    {
    if(d.allSame()) return std::sqrt(std::norm(d.val))*std::sqrt(minM(N.is));
    return vec_norm(d.store);
    }

Real
doTask(const NormNoScale& N, const ITCombiner& d) { return 0; }

//template<>
//void ITensor::
//scaleOutNorm()
//    {
//    auto nrm = doTask<Real>(NormNoScale{is_},store_);
//    //If norm already 1 return so
//    //we don't have to call MultReal
//    if(fabs(nrm-1.) < 1E-12) return;
//    if(nrm == 0)
//        {
//        scale_ = LogNumber(1.);
//        return;
//        }
//    doTask(MultReal{1./nrm},store_);
//    scale_ *= nrm;
//    }

//template<>
//void ITensor::
//equalizeScales(ITensor& other)
//    {
//    if(scale_.sign() != 0)
//        {
//        other.scaleTo(scale_);
//        }
//    else //*this is equivalent to zero
//        {
//        fill(0);
//        scale_ = other.scale_;
//        }
//    }

void
doTask(Conj, ITCplx& d) 
    { 
    auto i = MAKE_SAFE_PTR(d.istart(),d.csize());
    auto* ie = d.iend();
    for(; i != ie; ++i) *i *= -1;
    }

void
doTask(Conj, ITDiag<Cplx>& d) 
    { 
    if(d.allSame()) 
        {
        d.val = std::conj(d.val);
        }
    else
        {
        for(auto& el : d.store) 
            el = std::conj(el);
        }
    }

void
doTask(Conj,const ITReal& d) { }

void
doTask(Conj,const ITCombiner& d) { }

void
doTask(Conj,const ITDiag<Real>& d) { }

template<>
ITensor& ITensor::
conj()
    {
    doTask(Conj{},store_);
    return *this;
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
doTask(TakeReal,const ITDiag<Cplx>& d, ManagePtr& mp) 
    { 
    if(d.allSame()) mp.makeNewData<ITDiag<Real>>(d.val.real());
    else            
        {
        auto nd = mp.makeNewData<ITDiag<Real>>(d.size(),0.);
        for(auto i : index(d.store)) nd->store[i] = d.store[i].real();
        }
    }

void
doTask(TakeReal, const ITReal& ) { }

void
doTask(TakeReal, const ITDiag<Real>& ) { }


template<>
ITensor& ITensor::
takeReal()
    {
    doTask(TakeReal{},store_);
    return *this;
    }

void
doTask(TakeImag, const ITReal& d, ManagePtr& mp) 
    { 
    mp.makeNewData<ITReal>(d.size(),0);
    }

void
doTask(TakeImag,const ITCplx& d, ManagePtr& mp) 
    { 
    mp.makeNewData<ITReal>(d.istart(),d.iend());
    }

void
doTask(TakeImag,const ITDiag<Cplx>& d, ManagePtr& mp) 
    { 
    if(d.allSame()) mp.makeNewData<ITDiag<Real>>(d.val.imag());
    else            
        {
        auto nd = mp.makeNewData<ITDiag<Real>>(d.size(),0.);
        for(auto i : index(d.store)) nd->store[i] = d.store[i].imag();
        }
    }

template<typename T>
void
doTask(TakeImag,const T& d) { }

template<>
ITensor& ITensor::
takeImag()
    {
    doTask(TakeImag{},store_);
    return *this;
    }

ostream& 
operator<<(ostream & s, const ITensor& t)
    {
    s << "ITensor r=" << t.r() << ": " << t.inds() << "\n";
    if(!t) 
        {
        s << "{Storage is default constructed}\n";
        }
    else
        {
        //Checking whether std::ios::floatfield is set enables 
        //printing the contents of an ITensor when using the printf
        //format string %f (or another float-related format string)
        bool ff_set = (std::ios::floatfield & s.flags()) != 0;
        bool print_data = (ff_set || Global::printdat());
        doTask(PrintIT{s,t.scale(),t.inds(),print_data},t.cstore());
        }
    return s;
    }

void
doTask(PrintIT& P, const ITCombiner& d)
    {
    P.printInfo(d,"Combiner");
    }

void
doTask(PrintIT& P, const ITReal& d)
    {
    P.printInfo(d,"Dense Real",doTask(NormNoScale(P.is),d));
     
    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        detail::printVal(P.s,P.scalefac*d.store.front());
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

            detail::printVal(P.s,val);
            }
        }
    }

void
doTask(PrintIT& P, const ITCplx& d)
    {
    P.printInfo(d,"Dense Cplx",doTask(NormNoScale(P.is),d));

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        detail::printVal(P.s,P.scalefac*d.get(0));
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

            detail::printVal(P.s,val);
            }
        }
    }

template<typename T>
void
doTask(PrintIT& P, const ITDiag<T>& d)
    {
    constexpr auto type = std::is_same<T,Real>::value ? "Real" : "Cplx";
    P.printInfo(d,format("Diag %s%s",type,d.allSame()?", all same":""),
              doTask(NormNoScale(P.is),d));

    if(P.is.r() == 0) 
        {
        P.s << "  ";
        detail::printVal(P.s,P.scalefac*(d.empty() ? d.val : d.store.front()));
        return;
        }

    if(!P.print_data) return;

    auto size = minM(P.is);
    for(auto i : count(size))
        {
        auto val = P.scalefac*(d.allSame() ? d.val : d.store[i]);
        if(std::norm(val) > Global::printScale())
            {
            P.s << "(";
            for(size_t j = 1; j < P.is.size(); ++j)
                {
                P.s << (1+i) << ",";
                }
            P.s << (1+i) << ") ";
            detail::printVal(P.s,val);
            }
        }
    }

Cplx
quickranCplx() { return Cplx(detail::quickran(),detail::quickran()); }

ITensor
randomize(ITensor T, const Args& args)
    {
    if(args.getBool("Complex",false)) T.generate(quickranCplx);
    else                              T.generate(detail::quickran);
    return T;
    }

ITensor
matrixTensor(Mat&& M, const Index& i1, const Index& i2)
    {
    auto res = ITensor({i1,i2},ITReal{std::move(M.store())});
    M.clear();
    return res;
    }


Real
norm(const ITensor& T)
    {
#ifdef DEBUG
    if(!T) Error("ITensor is default initialized");
#endif
    return fabs(T.scale().real0()) *
           doTask<Real>(NormNoScale{T.inds()},T.cstore());
    }

ITensor
conj(ITensor T)
    {
    T.conj();
    return T;
    }

bool
doTask(CheckComplex,const ITCplx& d) { return true; }

bool
doTask(CheckComplex,const ITDiag<Cplx>& d) { return true; }

template<typename T>
bool
doTask(CheckComplex, const T& d) { return false; }

bool
isComplex(const ITensor& t)
    {
    return doTask<bool>(CheckComplex{},t.cstore());
    }

Cplx
doTask(SumEls, const ITReal& d) 
    { 
    Real sum = 0;
    for(const auto& elt : d)
        sum += elt;
    return sum;
    }

Cplx
doTask(SumEls, const ITCplx& d) 
    { 
    Real rsum = 0,
         isum = 0;
    auto p = MAKE_SAFE_PTR(d.rstart(),d.csize());
    for(; p != d.istart(); ++p) rsum += *p;
    p = MAKE_SAFE_PTR(d.istart(),d.csize());
    for(; p != d.iend(); ++p)   isum += *p;
    return Cplx(rsum,isum);
    }

template <class T>
Cplx
doTask(SumEls S, const ITDiag<T>& d) 
    { 
    if(d.allSame()) return Real(minM(S.is))*d.val;
    T sum = 0;
    for(const auto& elt : d.store)
        sum += elt;
    return sum;
    }

Cplx
sumelsC(const ITensor& t)
    {
    auto z = doTask<Cplx>(SumEls{t.inds()},t.cstore());
    return t.scale().real0()*z;
    }

Real
sumels(const ITensor& t)
    {
    auto z = sumelsC(t);
    if(z.imag() != 0) Error("ITensor has non-zero imaginary part, use sumelsC");
    return z.real();
    }

ITensor
combiner(std::vector<Index> inds, const Args& args)
    {
    if(inds.empty()) Error("No indices passed to combiner");
    long rm = 1;
    for(const auto& i : inds)
        {
        rm *= i.m();
        }
    //increase size by 1
    inds.push_back(Index());
    //shuffle contents to the end
    for(size_t j = inds.size()-1; j > 0; --j)
        {
        inds[j] = inds[j-1];
        }
    //create combined index
    auto cname = args.getString("IndexName","cmb");
    auto itype = getIndexType(args,"IndexType",Link);
    inds.front() = Index(cname,rm,itype);
    return ITensor(IndexSet(std::move(inds)),ITCombiner());
    }

ITensor
deltaTensor(const Index& i1, const Index& i2)
    {
#ifdef DEBUG
    if(i1.m() != i2.m()) Error("delta: indices must have same dimension");
#endif
    return ITensor({i1,i2},ITCombiner());
    }

enum class ITStorage { Null=0, Real=1, Cplx=2, Combiner=3, DiagReal=4, DiagCplx=5 }; 

template<typename T, typename... CtrArgs>
ITensor::storage_ptr
readType(std::istream& s, CtrArgs&&... args)
    {
    auto p = std::make_shared<ITDataType<T>>(std::forward<CtrArgs>(args)...);
    read(s,p->d);
    return p;
    }

void
read(std::istream& s, ITensor& T)
    {
    IndexSet is;
    read(s,is);
    LogNumber scale;
    read(s,scale);
    auto type = ITStorage::Null;
    s.read((char*)&type,sizeof(type));
    ITensor::storage_ptr p;
    if(type==ITStorage::Null) { /*intentionally left blank*/  }
    else if(type==ITStorage::Real) { p = readType<ITReal>(s); }
    else if(type==ITStorage::Cplx) { p = readType<ITCplx>(s); }
    else if(type==ITStorage::Combiner) { p = readType<ITCombiner>(s); }
    else if(type==ITStorage::DiagReal) { p = readType<ITDiag<Real>>(s); }
    else if(type==ITStorage::DiagCplx) { p = readType<ITDiag<Cplx>>(s); }
    else
        {
        Error("Unrecognized type when reading ITensor from istream");
        }
    T = ITensor(std::move(is),std::move(p),scale);
    }

template<class T>
void
writeType(std::ostream& s, ITStorage type, const T& data)
    {
    s.write((char*)&type,sizeof(type));
    write(s,data); 
    }

void
doTask(Write& W, const ITReal& d) 
    { 
    writeType(W.s,ITStorage::Real,d); 
    }

void
doTask(Write& W, const ITCplx& d) 
    { 
    writeType(W.s,ITStorage::Cplx,d); 
    }

void
doTask(Write& W, const ITCombiner& d) 
    { 
    writeType(W.s,ITStorage::Combiner,d); 
    }

void
doTask(Write& W, const ITDiag<Real>& d)
    { 
    writeType(W.s,ITStorage::DiagReal,d); 
    }

void
doTask(Write& W, const ITDiag<Cplx>& d)
    { 
    writeType(W.s,ITStorage::DiagCplx,d); 
    }

void
write(std::ostream& s, const ITensor& T)
    {
    write(s,T.inds());
    write(s,T.scale());
    if(T) 
        {
        doTask(Write{s},T.cstore());
        }
    else 
        {
        auto type = ITStorage::Null;
        s.write((char*)&type,sizeof(type));
        }
    }

} //namespace itensor
