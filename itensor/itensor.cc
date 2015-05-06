//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor.h"
#include "lapack_wrap.h"
#include "detail/printing.h"
#include "detail/gcounter.h"
#include "contract.h"
#include "count.h"

using std::array;
using std::ostream;
using std::vector;
using std::shared_ptr;
using std::make_shared;
using std::cout;
using std::endl;

namespace itensor {

//
// ITensor Constructors
//


ITensor::
ITensor(const Index& i1) 
    :
    is_(i1),
    store_(std::make_shared<ITReal>(i1.m(),0.)),
    scale_(1.)
	{ }


ITensor::
ITensor(const Index& i1,const Index& i2) 
    :
    is_(i1,i2),
    store_(std::make_shared<ITReal>(i1.m()*i2.m(),0.)),
    scale_(1.)
	{ }
    
ITensor::
ITensor(Complex val) 
    :
    scale_(1.)
    { 
    //if(val.imag() == 0)
    //    store_ = std::make_shared<ITReal>(1,val.real());
    //else
    //    store_ = std::make_shared<ITCplx>(1,val);
    if(val.imag() == 0)
        store_ = std::make_shared<ITDiag<Real>>(val.real());
    else
        store_ = std::make_shared<ITDiag<Complex>>(val);
    }

//ITensor::
//ITensor(IndexSet&& iset,
//        storage_ptr&& data,
//        const LogNumber& scale = 1)
//    :
//    is_(std::move(iset)),
//    store_(std::move(data)),
//    scale_(scale)
//    { }

ITensor::
ITensor(const IndexSet& is)
    :
    is_(is),
    store_(std::make_shared<ITReal>(area(is_),0.)),
    scale_(1.)
	{ }

struct CopyElems : public RegisterFunc<CopyElems>
    {
    CopyElems() { }

    void
    operator()(ITReal& d1,
               const ITReal& d2)
        {
        std::copy(d2.begin(),d2.end(),d1.begin());
        }

    template<typename T>
    void
    operator()(ITDiag<T>& d1,
               const ITDiag<T>& d2)
        {
        std::copy(d2.begin(),d2.end(),d1.begin());
        }
    };

//ITensor::
//ITensor(const Index& i1,
//        const Index& i2,
//        const MatrixRef& M)
//    :
//    is_(i1,i2),
//    scale_(1.)
//    {
//	if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
//	    Error("Mismatch of Index sizes and matrix.");
//
//    if(i1 == is_[0])
//        {
//        Matrix Mt = M.t();
//        VectorRef vref = Mt.TreatAsVector(); 
//        store_ = make_newdata<ITReal>(vref.begin(),vref.end());
//        }
//    else
//        {
//        VectorRef vref = M.TreatAsVector(); 
//        store_ = make_newdata<ITReal>(vref.begin(),vref.end());
//        }
//    }

ITensor::
ITensor(IndexSet iset,
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

struct Contract : RegisterFunc<Contract>
    {
    private:
    const Label &Lind_,
                &Rind_;

    const IndexSet &Lis_,
                   &Ris_;

    //New IndexSet
    IndexSet Nis_;
    Real scalefac_ = 0;

    public:

    Contract(const IndexSet& Lis,
             const Label& Lind,
             const IndexSet& Ris,
             const Label& Rind)
        :
        Lind_(Lind),
        Rind_(Rind),
        Lis_(Lis),
        Ris_(Ris)
        { }

    IndexSet
    newIndexSet() { return std::move(Nis_); }
    Real
    scalefac() { return scalefac_; }

    void
    operator()(const ITReal& a1,
               const ITReal& a2);

    void
    operator()(const ITCplx& a1,
               const ITCplx& a2);


    void
    operator()(const ITReal& a1,
               const ITCplx& a2)
        {
        realCplx(a1,Lis_,Lind_,a2,Ris_,Rind_);
        }

    void
    operator()(const ITCplx& a1,
               const ITReal& a2)
        {
        realCplx(a2,Ris_,Rind_,a1,Lis_,Lind_);
        }

    void
    operator()(const ITReal& d,
               const ITCombiner& C)
        {
        combine(d,Lis_,Ris_);
        }
    void
    operator()(const ITCombiner& C,
               const ITReal& d)
        { 
        combine(d,Ris_,Lis_);
        if(!newData()) assignPointerRtoL();
        }

    void
    operator()(const ITDiag<Real>& d,
               const ITReal& t)
        {
        diagDense(d,Lis_,Lind_,t,Ris_,Rind_);
        }
    void
    operator()(const ITReal& t,
               const ITDiag<Real>& d)
        { 
        diagDense(d,Ris_,Rind_,t,Lis_,Lind_);
        }

    private:

    void
    combine(const ITReal& d,
            const IndexSet& dis,
            const IndexSet& Cis);

    void
    diagDense(const ITDiag<Real>& d,
              const IndexSet& dis,
              const Label& dind,
              const ITReal& t,
              const IndexSet& tis,
              const Label& tind);

    void
    realCplx(const ITReal& R,
             const IndexSet& ris,
             const Label& rind,
             const ITCplx& C,
             const IndexSet& cis,
             const Label& cind);

    template<typename Data>
    void
    computeScalefac(Data& dat)
        {
        scalefac_ = 0;
        for(auto elt : dat) scalefac_ += elt*elt;
        scalefac_ = std::sqrt(scalefac_);
        if(scalefac_ == 0) return;
        for(auto& elt : dat) elt /= scalefac_;
        }
 
    };

void Contract::
operator()(const ITCplx& a1,
           const ITCplx& a2)
    {
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    const bool sortInds = false; //whether to sort indices of result
    contractIS(Lis_,Lind_,Ris_,Rind_,Nis_,sortInds);
    
    Label Nind(Nis_.r(),0);
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(Lis_,Nis_[i]);
        if(j >= 0)
            {
            Nind[i] = Lind_[j];
            }
        else
            {
            j = findindex(Ris_,Nis_[i]);
            Nind[i] = Rind_[j];
            }
        }

    auto rsize = area(Nis_);
    auto nd = makeNewData<ITCplx>(rsize,0.);

    auto t1r = make_tensorref(a1.rstart(),Lis_),
         t1i = make_tensorref(a1.istart(),Lis_),
         t2r = make_tensorref(a2.rstart(),Ris_),
         t2i = make_tensorref(a2.istart(),Ris_);
    auto trr = make_tensorref(nd->rstart(),Nis_),
         tri = make_tensorref(nd->istart(),Nis_);

    contractloop(t1i,Lind_,t2i,Rind_,trr,Nind);
    for(auto p = nd->rstart(); p < nd->istart(); ++p) *p *= -1;
    contractloop(t1r,Lind_,t2r,Rind_,trr,Nind);

    contractloop(t1i,Lind_,t2r,Rind_,tri,Nind);
    contractloop(t1r,Lind_,t2i,Rind_,tri,Nind);

    if(rsize > 1) computeScalefac(*nd);
    }

void Contract::
realCplx(const ITReal& R,
         const IndexSet& ris,
         const Label& rind,
         const ITCplx& C,
         const IndexSet& cis,
         const Label& cind)
    {
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    const bool sortInds = false; //whether to sort indices of result
    contractIS(ris,rind,cis,cind,Nis_,sortInds);
    
    Label Nind(Nis_.r(),0);
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(ris,Nis_[i]);
        if(j >= 0)
            {
            Nind[i] = rind[j];
            }
        else
            {
            j = findindex(cis,Nis_[i]);
            Nind[i] = cind[j];
            }
        }

    auto rsize = area(Nis_);
    auto nd = makeNewData<ITCplx>(rsize,0.);

    auto t1 = make_tensorref(R.data(),ris),
         t2r = make_tensorref(C.rstart(),cis),
         t2i = make_tensorref(C.istart(),cis);
    auto trr = make_tensorref(nd->rstart(),Nis_),
         tri = make_tensorref(nd->istart(),Nis_);

    contractloop(t1,rind,t2r,cind,trr,Nind);
    contractloop(t1,rind,t2i,cind,tri,Nind);

    if(rsize > 1) computeScalefac(*nd);
    }

void Contract::
operator()(const ITReal& a1,
           const ITReal& a2)
    {
    //Optimization TODO:
    //  Test different scenarios where having sortInds=true or false
    //  can improve performance. Having sorted inds can make adding
    //  quicker and let contractloop run in parallel more often in principle.
    const bool sortInds = false; //whether to sort indices of result
    contractIS(Lis_,Lind_,Ris_,Rind_,Nis_,sortInds);
    
    Label Nind(Nis_.r(),0);
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(Lis_,Nis_[i]);
        if(j >= 0)
            {
            Nind[i] = Lind_[j];
            }
        else
            {
            j = findindex(Ris_,Nis_[i]);
            Nind[i] = Rind_[j];
            }
        }

    //PRI(Lind_);
    //PRI(Rind_);
    //PRI(Nind);

    auto rsize = area(Nis_);
    auto nd = makeNewData<ITReal>(rsize,0.);
    auto t1 = make_tensorref(a1.data(),Lis_),
         t2 = make_tensorref(a2.data(),Ris_);
    auto tr = make_tensorref(nd->data(),Nis_);
    contractloop(t1,Lind_,t2,Rind_,tr,Nind);

    if(rsize > 1) computeScalefac(*nd);
    }

void Contract::
diagDense(const ITDiag<Real>& d,
          const IndexSet& dis,
          const Label& dind,
          const ITReal& t,
          const IndexSet& tis,
          const Label& tind)
    {
    contractIS(Lis_,Lind_,Ris_,Rind_,Nis_,true);

    long t_cstride = 0; //total t-stride of contracted inds of t
    size_t ntu = 0; //number uncontracted inds of t
    assert(int(tind.size()) == tis.size());
    for(size_t j = 0; j < tind.size(); ++j)
        {
        //if index j is contracted, add its stride to t_cstride:
        if(tind[j] < 0) t_cstride += tis.stride(j);
        else            ++ntu;
        }

    long d_ustride = 0; //total result-stride of uncontracted inds of d
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(dis,Nis_[i]);
        if(j >= 0) d_ustride += Nis_.stride(i);
        }

    auto dsize = size_t(minM(dis));

    if(ntu > 0)
        {
        vector<long> tstride(ntu,0),
                     rstride(ntu,0);
        detail::GCounter C(0,ntu,0);
        size_t n = 0;
        for(size_t j = 0; j < tind.size(); ++j)
            {
            if(tind[j] > 0)
                {
#ifdef DEBUG
                if(n >= ntu) Error("n out of range");
#endif
                C.setInd(n,0,tis.dim(j)-1);
                tstride.at(n) = tis.stride(j);
                auto k = findindex(Nis_,tis[j]);
#ifdef DEBUG
                if(k < 0) Error("Index not found");
#endif
                rstride.at(n) = Nis_.stride(k);
                ++n;
                }
            }
        auto nd = makeNewData<ITReal>(area(Nis_),0.);
        auto *pr = nd->data();
        const auto *pt = t.data();

        if(d.allSame())
            {
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(size_t i = 0; i < ntu; ++i)
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(long J = 0; J < dsize; ++J)
                    {
                    pr[J*d_ustride+roffset] += d.val*pt[J*t_cstride+toffset];
                    }
                }
            }
        else
            {
            auto* pd = d.data();
            assert(d.size() == dsize);
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(size_t i = 0; i < ntu; ++i)
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(size_t J = 0; J < dsize; ++J)
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
            const auto *pt = t.data();
            if(d.allSame())
                {
                for(size_t J = 0; J < dsize; ++J)
                    val += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(dsize == d.size());
                auto *pd = d.data();
                for(size_t J = 0; J < dsize; ++J)
                    val += pd[J]*pt[J*t_cstride];
                }
            makeNewData<ITDiag<Real>>(val);
            }
        else //some of d's inds uncontracted
            {
            // o element-wise product of d's data and t's diagonal
            auto nd = makeNewData<ITDiag<Real>>(dsize,0.);
            auto *pr = nd->data();
            const auto *pt = t.data();
            if(d.allSame())
                {
                for(size_t J = 0; J < dsize; ++J)
                    pr[J] += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(dsize == d.size());
                auto *pd = d.data();
                for(size_t J = 0; J < dsize; ++J)
                    pr[J] += pd[J]*pt[J*t_cstride];
                }
            }
        }
    }

void Contract::
combine(const ITReal& d,
        const IndexSet& dis,
        const IndexSet& Cis)
    {
    //TODO: try to make use of Lind,Rind label vectors
    //      to simplify combine logic
    const auto& cind = Cis[0];
    int jc = findindex(dis,cind);
    if(jc >= 0) //has cind
        {
        //dis has cind, replace with other inds
        vector<Index> newind;
        newind.reserve(dis.r()+Cis.r()-2);
        for(int j = 0; j < dis.r(); ++j)
            if(j == jc)
                {
                for(int k = 1; k < Cis.size(); ++k)
                    newind.push_back(Cis[k]);
                }
            else
                {
                newind.push_back(dis[j]);
                }
        Nis_ = IndexSet(move(newind));
        }
    else
        {
        //dis doesn't have cind, replace
        //Cis[1], Cis[2], ... with cind
        //may need to permute
        int J1 = findindex(dis,Cis[1]);
        if(J1 < 0) 
            {
            println("IndexSet of dense tensor = \n",dis);
            println("IndexSet of combiner/delta = \n",Cis);
            Error("No contracted indices in combiner-tensor product");
            }
        //Check if Cis[1],Cis[2],... are grouped together (contiguous)
        bool contig = true;
        for(int j = J1+1, c = 2; c < Cis.r() && j < dis.r(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                contig = false;
                break;
                }
        if(contig)
            {
            vector<Index> newind;
            newind.reserve(dis.r()-Cis.r()+2);
            for(int j = 0; j < J1; ++j) 
                newind.push_back(dis[j]);
            newind.push_back(cind);
            for(int j = J1+Cis.r()-1; j < dis.r(); ++j) 
                newind.push_back(dis[j]);
            assert(newind.size() == dis.r()-Cis.r()+2);
            Nis_ = IndexSet(move(newind));
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
            for(int c = 1; c < Cis.r(); ++c)
                {
                int j = findindex(dis,Cis[c]);
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
            vector<Index> newind;
            vector<long> pdims(dis.r(),-1);
            newind.reserve(dis.r()-Cis.r()+1);
            newind.push_back(cind);
            for(int j = 0; j < dis.r(); ++j)
                {
                if(P.dest(j) == -1) 
                    {
                    P.setFromTo(j,ni++);
                    newind.push_back(dis[j]);
                    }
                pdims[j] = dis[P.dest(j)].m();
                }
            Range rr(pdims);
            Nis_ = IndexSet(move(newind));
            auto nd = makeNewData<ITReal>(area(Nis_));
            auto td = make_tensorref(d.data(),dis);
            auto tr = make_tensorref(nd->data(),rr);
            permute(td,P,tr);
            }
        }
    }


ITensor& ITensor::
operator*=(const ITensor& other)
    {
    if(!(*this) || !other)
        Error("Default constructed ITensor in product");

    if(this == &other)
        return operator=( ITensor(sqr(norm(*this))) );

    const auto& Lis = is_;
    const auto& Ris = other.is_;

    Label Lind,
          Rind;
    auto ncont = computeLabels(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind);
    auto nuniq = Lis.r()+Ris.r()-2*ncont;

    //Check if other is a scalar (modulo m==1 inds)
    if(Ris.rn() == 0)
        {
        operator*=(other.cplx());
        is_ = IndexSet(computeNewInds(Lis,Lind,Ris,Rind,nuniq));
        return *this;
        }
    //Check if this is a scalar (modulo m==1 inds)
    if(Lis.rn() == 0)
        {
        auto newind = computeNewInds(Lis,Lind,Ris,Rind,nuniq);
        operator=(other*cplx());
        is_ = IndexSet(std::move(newind));
        return *this;
        }

    auto C = applyFunc<Contract>(store_,other.store_,Lis,Lind,Ris,Rind);

    is_ = C.newIndexSet();

    scale_ *= other.scale_;
    if(C.scalefac() > 0) scale_ *= C.scalefac();

    return *this;
    }

ITensor& ITensor::
operator*=(Real fac)
    {
    if(fac == 0)
        {
        fill(0);
        return *this;
        }
    scale_ *= fac;
    return *this;
    }

struct MultComplex : RegisterFunc<MultComplex>
    {
    private:
    Complex z_;
    public:
    MultComplex(Complex z) : z_(z) { }

    void
    operator()(const ITReal& d)
        {
        auto nd = makeNewData<ITCplx>(d);
        *nd *= z_;
        }
    void
    operator()(ITCplx& d) { d *= z_; }
    };

ITensor& ITensor::
operator*=(Complex z)
    {
    if(z.imag() == 0) return operator*=(z.real());
    applyFunc<MultComplex>(store_,z);
    return *this;
    }

class PlusEQ : public RegisterFunc<PlusEQ>
    {
    Real fac_;
    const Permutation *P_ = nullptr;
    const IndexSet *is1_ = nullptr,
                   *is2_ = nullptr;
    public:
    using permutation = Permutation;

    PlusEQ(Real fac)
        :
        fac_(fac)
        { }

    PlusEQ(const Permutation& P,
           const IndexSet& is1,
           const IndexSet& is2,
           Real fac)
        :
        fac_(fac),
        P_(&P),
        is1_(&is1),
        is2_(&is2)
        { }

    void
    operator()(ITReal& a1,
               const ITReal& a2);

    void
    operator()(ITDiag<Real>& a1,
               const ITDiag<Real>& a2);

    void
    operator()(ITReal& a1,
               const ITCplx& a2)
        {
        Error("Real + Complex not implemented");
        //auto np = make_newdata<ITCplx>(a1);
        //operator()(*np,a2);
        }

    void
    operator()(ITCplx& a1,
               const ITReal& a2)
        {
        Error("Complex + Real not implemented");
        //ITCplx a2c(a2);
        //operator()(a1,a2c);
        }
    };

void
plusEqData(Real fac, Real *d1, const Real *d2, LAPACK_INT size)
    {
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,d2,&inc,d1,&inc);
    }

void PlusEQ::
operator()(ITReal& a1,
           const ITReal& a2)
    {
#ifdef DEBUG
    if(a1.size() != a2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(!P_)
        {
        plusEqData(fac_,a1.data(),a2.data(),a1.size());
        }
    else
        {
        auto ref1 = tensorref<Real,IndexSet>(a1.data(),*is1_),
             ref2 = tensorref<Real,IndexSet>(a2.data(),*is2_);
        auto f = fac_;
        auto add = [f](Real& r1, Real r2) { r1 += f*r2; };
        permute(ref2,*P_,ref1,add);
        }
    }

void PlusEQ::
operator()(ITDiag<Real>& a1,
           const ITDiag<Real>& a2)
    {
#ifdef DEBUG
    if(a1.size() != a2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(a1.allSame() || a2.allSame()) Error("ITDiag plusEq allSame case not implemented");
    plusEqData(fac_,a1.data(),a2.data(),a1.size());
    }

ITensor& ITensor::
operator+=(const ITensor& other)
    {
    if(!*this) Error("Calling += on default constructed ITensor");
    if(!other) Error("Right-hand-side of ITensor += is default constructed");
    if(this == &other) return operator*=(2.);
    if(this->scale_.isZero()) return operator=(other);

    PlusEQ::permutation P(is_.size());
#ifdef DEBUG
    try {
        detail::calc_permutation(other.is_,is_,P);
        }
    catch(const std::exception& e)
        {
        Print(*this);
        Print(other);
        Error("ITensor::operator+=: different Index structure");
        }
#else
    detail::calc_permutation(other.is_,is_,P);
#endif


    Real scalefac = 1;
    if(scale_.magnitudeLessThan(other.scale_)) 
        {
        this->scaleTo(other.scale_); 
        }
    else
        {
        scalefac = (other.scale_/scale_).real();
        }

    if(isTrivial(P))
        {
        applyFunc<PlusEQ>(store_,other.store_,scalefac);
        }
    else
        {
        applyFunc<PlusEQ>(store_,other.store_,P,is_,other.is_,scalefac);
        }

    return *this;
    } 

ITensor& ITensor::
operator-=(const ITensor& other)
    {
    if(this == &other) 
        { 
        scale_ = 0; 
        fill(0);
        return *this; 
        }
    scale_.negate();
    operator+=(other); 
    scale_.negate();
    return *this; 
    }

class FillReal : public RegisterFunc<FillReal>
    {
    Real r_;
    public:
    FillReal(Real r) : r_(r) { }
    void
    operator()(ITReal& d) const
        {
        std::fill(d.begin(),d.end(),r_);
        }
    void
    operator()(const ITCplx& d)
        {
        makeNewData<ITReal>(d.csize(),r_);
        }
    template<typename T>
    void
    operator()(const ITDiag<T>& d)
        {
        makeNewData<ITDiag<Real>>(r_);
        }
    };

class FillCplx : public RegisterFunc<FillCplx>
    {
    Complex z_;
    public:
    FillCplx(Complex z) : z_(z) { }
    void
    operator()(const ITReal& d)
        {
        makeNewData<ITCplx>(d.size(),z_);
        }
    void
    operator()(ITCplx& d) const
        {
        d.fill(z_);
        }
    template<typename T>
    void
    operator()(const ITDiag<T>& d)
        {
        makeNewData<ITDiag<Complex>>(z_);
        }
    };

ITensor& ITensor::
fill(Complex z)
    {
    if(!(*this)) return *this;
    scale_ = LogNumber(1.);
    if(z.imag() == 0)
        applyFunc<FillReal>(store_,z.real());
    else
        applyFunc<FillCplx>(store_,z);
    return *this;
    }

class MultReal : public RegisterFunc<MultReal>
    {
    Real r_;
    public:
    MultReal(Real r)
        : r_(r)
        { }

    void
    operator()(ITReal& d) const
        {
        //use BLAS algorithm?
        for(auto& elt : d) elt *= r_;
        }

    void
    operator()(ITCplx& d) const
        {
        for(auto& elt : d) elt *= r_;
        }

    template<typename T>
    void
    operator()(ITDiag<T>& d) const
        {
        d.val *= r_;
        //use BLAS algorithm?
        for(auto& elt : d.store) elt *= r_;
        }
    };

void ITensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    applyFunc<MultReal>(store_,scale_.real0());
    scale_ = newscale;
    }


class NormNoScale : public RegisterFunc<NormNoScale,Real>
    {
    const IndexSet& is_;
    public:

    NormNoScale(const IndexSet& is) : is_(is) { }

    Real
    operator()(const ITReal& d) { return vec_norm(d); }

    Real
    operator()(const ITCplx& d) { return vec_norm(d); }

    template<typename T>
    Real
    operator()(const ITDiag<T>& d)
        {
        if(d.allSame()) return std::sqrt(std::norm(d.val))*std::sqrt(minM(is_));
        return vec_norm(d.store);
        }

    private:

    template<typename Container>
    Real
    vec_norm(const Container& v)
        {
        Real nrm = 0;
        for(const auto& elt : v)
            nrm += std::norm(elt); //conj(elt)*elt
        return std::sqrt(nrm);
        }
    };

void ITensor::
scaleOutNorm()
    {
    auto nrm = applyFunc<NormNoScale>(store_,is_);
    //If norm already 1 return so
    //we don't have to call MultReal
    if(fabs(nrm-1.) < 1E-12) return;
    if(nrm == 0)
        {
        scale_ = LogNumber(1.);
        return;
        }
    applyFunc<MultReal>(store_,1./nrm);
    scale_ *= nrm;
    }

void ITensor::
equalizeScales(ITensor& other)
    {
    if(scale_.sign() != 0)
        {
        other.scaleTo(scale_);
        }
    else //*this is equivalent to zero
        {
        fill(0);
        scale_ = other.scale_;
        }
    }

struct Conj : RegisterFunc<Conj>
    {
    void
    operator()(const ITCplx& cd) 
        { 
        auto& d = modifyData(cd);
        auto* i = d.istart();
        auto* ie = i+d.csize();
        for(; i < ie; ++i) *i *= -1;
        }
    void
    operator()(const ITDiag<Complex>& cd) 
        { 
        auto& d = modifyData(cd);
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
    template<typename T>
    void
    operator()(const T& d) { }
    };

ITensor& ITensor::
conj()
    {
    applyFunc<Conj>(store_);
    return *this;
    }

struct TakeReal : RegisterFunc<TakeReal>
    {
    void
    operator()(const ITCplx& d) 
        { 
        makeNewData<ITReal>(d.rstart(),d.istart());
        }
    void
    operator()(const ITDiag<Complex>& d) 
        { 
        if(d.allSame()) makeNewData<ITDiag<Real>>(d.val.real());
        else            
            {
            auto nd = makeNewData<ITDiag<Real>>(d.size(),0.);
            for(auto i : index(d.store)) nd->store[i] = d.store[i].real();
            }
        }
    template<typename T>
    void
    operator()(const T& d) { }
    };
ITensor& ITensor::
takeReal()
    {
    applyFunc<TakeReal>(store_);
    return *this;
    }

struct TakeImag : RegisterFunc<TakeImag>
    {
    void
    operator()(const ITCplx& d) 
        { 
        makeNewData<ITReal>(d.istart(),d.iend());
        }
    void
    operator()(const ITDiag<Complex>& d) 
        { 
        if(d.allSame()) makeNewData<ITDiag<Real>>(d.val.imag());
        else            
            {
            auto nd = makeNewData<ITDiag<Real>>(d.size(),0.);
            for(auto i : index(d.store)) nd->store[i] = d.store[i].imag();
            }
        }
    template<typename T>
    void
    operator()(const T& d) { }
    };
ITensor& ITensor::
takeImag()
    {
    applyFunc<TakeImag>(store_);
    return *this;
    }

struct PrintIT : RegisterFunc<PrintIT>
    {
    std::ostream& s_;
    const LogNumber& x_;
    const IndexSet& is_;

    PrintIT(std::ostream& s,
            const LogNumber& x,
            const IndexSet& is)
        : s_(s), x_(x), is_(is)
        { }

    void
    operator()(const ITReal& d) const;

    void
    operator()(const ITCplx& d) const;

    template<typename T>
    void
    operator()(const ITDiag<T>& d) const;

    void
    operator()(const ITCombiner& d) const { s_ << " Combiner}\n"; }
    };

void PrintIT::
operator()(const ITReal& d) const
    {
    s_ << " (Dense Real)}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) 
        {
        s_ << "  ";
        detail::printVal(s_,scalefac*d.store.front());
        return;
        }

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,is_.dim(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = scalefac*d[ind(is_,gc.i)];
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            detail::printVal(s_,val);
            }
        }
    }

void PrintIT::
operator()(const ITCplx& d) const
    {
    s_ << " (Dense Cplx)}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) 
        {
        s_ << "  ";
        detail::printVal(s_,scalefac*d.get(0));
        return;
        }

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,is_.dim(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = scalefac*d.get(ind(is_,gc.i));
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            detail::printVal(s_,val);
            }
        }
    }

template<typename T>
void PrintIT::
operator()(const ITDiag<T>& d) const
    {
    constexpr auto type = std::is_same<T,Real>::value ? "Real" : "Cplx";
    auto allsame = d.allSame();
    s_ << " (Diag " << type << (allsame ? ",all same)" : ")") << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    if(is_.r() == 0) 
        {
        s_ << "  ";
        detail::printVal(s_,scalefac*(d.empty() ? d.val : d.store.front()));
        return;
        }

    auto size = minM(is_);
    for(size_t i = 0; i < size; ++i)
        {
        auto val = scalefac*(allsame ? d.val : d.store[i]);
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(int j = 1; j < is_.size(); ++j)
                {
                s_ << (1+i) << ",";
                }
            s_ << (1+i) << ") ";
            detail::printVal(s_,val);
            }
        }
    }
//template void PrintIT::operator()(const ITDiag<Real>& d) const;
//template void PrintIT::operator()(const ITDiag<Complex>& d) const;

ostream& 
operator<<(ostream & s, const ITensor& t)
    {
    s << "ITensor r=" << t.r() << ": ";
    s << t.inds() << "\n";
    s << "  {log(scale)[incl in elems]=" << t.scale().logNum();

    //Checking whether std::ios::floatfield is set enables 
    //printing the contents of an ITensor when using the printf
    //format string %f (or another float-related format string)
    const bool ff_set = (std::ios::floatfield & s.flags()) != 0;

    if(ff_set || Global::printdat())
        {
        if(t) applyFunc<PrintIT>(t.data(),s,t.scale(),t.inds());
        else           s << " (default constructed)}\n";
        }
    else
        {
        s << "}";
        }
    return s;
    }


Complex
quickranCplx() { return Complex(detail::quickran(),detail::quickran()); }

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
    return T.scale().real0() *
           applyFunc<NormNoScale>(T.data(),T.inds());
    }


ITensor
conj(ITensor T)
    {
    T.conj();
    return T;
    }


struct CheckComplex : RegisterFunc<CheckComplex,bool>
    {
    bool
    operator()(const ITCplx& d) { return true; }
    bool
    operator()(const ITDiag<Complex>& d) { return true; }

    template<typename T>
    bool
    operator()(const T& d) { return false; }
    };

bool
isComplex(const ITensor& t)
    {
    return applyFunc<CheckComplex>(t.data());
    }

class SumEls : public RegisterFunc<SumEls,Complex>
    {
    const IndexSet& is_;
    public:

    SumEls(const IndexSet& is) : is_(is) { }

    Complex
    operator()(const ITReal& d) 
        { 
        Real sum = 0;
        for(const auto& elt : d)
            sum += elt;
        return sum;
        }

    Complex
    operator()(const ITCplx& d) 
        { 
        Real rsum = 0,
             isum = 0;
        auto* p = d.rstart();
        for(; p < d.istart(); ++p) rsum += *p;
        for(; p < d.iend(); ++p)   isum += *p;
        return Complex(rsum,isum);
        }

    template <class T>
    Complex
    operator()(const ITDiag<T>& d) 
        { 
        if(d.allSame()) return Real(minM(is_))*d.val;
        T sum = 0;
        for(const auto& elt : d.store)
            sum += elt;
        return sum;
        }
    };

Complex
sumelsC(const ITensor& t)
    {
    auto z = Complex(applyFunc<SumEls>(t.data(),t.inds()));
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
combiner(std::vector<Index> inds)
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
    inds.front() = Index("cmb",rm);
    return ITensor(IndexSet(std::move(inds)),ITCombiner());
    }

ITensor
delta(const Index& i1, const Index& i2)
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
    auto p = std::make_shared<T>(std::forward<CtrArgs>(args)...);
    read(s,*p);
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
    else if(type==ITStorage::Real) { println("Reading ITensor of type ITReal"); p = readType<ITReal>(s); }
    else if(type==ITStorage::Cplx) { p = readType<ITCplx>(s); }
    else if(type==ITStorage::Combiner) { p = readType<ITCombiner>(s); }
    else if(type==ITStorage::DiagReal) { p = readType<ITDiag<Real>>(s); }
    else if(type==ITStorage::DiagCplx) { p = readType<ITDiag<Cplx>>(s); }
    else
        {
        Error("Unrecognized type when reading ITensor from istream");
        }
    //println("Data before returning from read:");
    //applyFunc<PrintIT>(p,std::cout,scale,is);
    T = ITensor(std::move(is),std::move(p),scale);
    }

template<class T>
void
writeType(std::ostream& s, ITStorage type, const T& data)
    {
    s.write((char*)&type,sizeof(type));
    write(s,data); 
    }

struct Write : RegisterFunc<Write>
    {
    std::ostream& s_;
    Write(std::ostream& s) : s_(s) { }

    void
    operator()(const ITReal& d)
        { 
        writeType(s_,ITStorage::Real,d);
        }

    void
    operator()(const ITCplx& d)
        { 
        writeType(s_,ITStorage::Cplx,d);
        }

    void
    operator()(const ITCombiner& d)
        { 
        writeType(s_,ITStorage::Combiner,d);
        }

    void
    operator()(const ITDiag<Real>& d)
        { 
        writeType(s_,ITStorage::DiagReal,d);
        }

    void
    operator()(const ITDiag<Cplx>& d)
        { 
        writeType(s_,ITStorage::DiagCplx,d);
        }
    };

void
write(std::ostream& s, const ITensor& T)
    {
    write(s,T.inds());
    write(s,T.scale());
    if(T) 
        {
        applyFunc<Write>(T.data(),s);
        }
    else 
        {
        auto type = ITStorage::Null;
        s.write((char*)&type,sizeof(type));
        }
    }

};
