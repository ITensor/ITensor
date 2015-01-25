//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor.h"

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
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(i1.m(),0.))
	{ }


ITensor::
ITensor(const Index& i1,const Index& i2) 
    :
    is_(i1,i2),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(i1.m()*i2.m(),0.))
	{ }
    
ITensor::
ITensor(Complex val) 
    :
    scale_(1.)
    { 
    if(val.imag() == 0)
        store_ = std::make_shared<ITDense<Real>>(1,val.real());
    else
        store_ = std::make_shared<ITDense<Complex>>(1,val);
    }

ITensor::
ITensor(IndexSet&& iset,
        NewData nd,
        LogNumber scale)
    :
    is_(iset),
    scale_(scale),
    store_(std::move(nd))
    { }

ITensor::
ITensor(const IndexSet& is)
    :
    is_(is),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(area(is_),0.))
	{ }

ITensor::
ITensor(const IndexSet& is,
        const VectorRef& v)
    :
    is_(is),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(v.begin(),v.end()))
	{ }

struct CopyElems
    {
    CopyElems() { }

    template<typename T>
    NewData
    operator()(ITDense<T>& d1,
               const ITDense<T>& d2)
        {
        std::copy(d2.data.begin(),d2.data.end(),d1.data.begin());
        return NewData();
        }

    template<typename T>
    NewData
    operator()(ITDiag<T>& d1,
               const ITDiag<T>& d2)
        {
        std::copy(d2.data.begin(),d2.data.end(),d1.data.begin());
        return NewData();
        }

    template<typename T1, typename T2>
    NewData
    operator()(T1& d1,
               const T2& d2) { Error("Not implemented"); return NewData(); }
    };

ITensor::
ITensor(const IndexSet& is,
        const ITensor& t)
    :
    is_(is),
    scale_(t.scale_),
    store_(std::make_shared<ITDense<Real>>(area(is_),0.))
    {
    Error("ITensor(IndexSet,ITensor) constructor currently broken due to automatic sorting of Indices by IndexSet");
    applyFunc<CopyElems>(store_,t.store_);
    }

ITensor::
ITensor(const Index& i1,
        const Index& i2,
        const MatrixRef& M)
    :
    is_(i1,i2),
    scale_(1.)
    {
	if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
	    Error("Mismatch of Index sizes and matrix.");

    if(i1 == is_[0])
        {
        Matrix Mt = M.t();
        VectorRef vref = Mt.TreatAsVector(); 
        store_ = make_newdata<ITDense<Real>>(vref.begin(),vref.end());
        }
    else
        {
        VectorRef vref = M.TreatAsVector(); 
        store_ = make_newdata<ITDense<Real>>(vref.begin(),vref.end());
        }
    }


//class Reshape
//    {
//    const Permutation& P_;
//    const IndexSet& is_;
//    public:
//    Reshape(const Permutation& P,
//            const IndexSet& is)
//        : P_(P), is_(is)
//        { }
//
//    NewData
//    operator()(ITDense<Real>& t1, 
//               const ITDense<Real>& t2)
//        {
//        auto v = std::vector<Real>();
//        reshape(P_,is_,t2.data(),v);
//        return NewData();
//        }
//
//    template<typename T1, typename T2>
//    NewData
//    operator()(T1& t1, const T2& t2)
//        {
//        Error("Reshape not implemented for ITData types");
//        return NewData();
//        }
//    };
//
//ITensor::
//ITensor(const IndexSet& is,
//        const ITensor& t,
//        const Permutation& P)
//    :
//    is_(is_),
//    scale_(t.scale_)
//    {
//    if(isTrivial(P)) 
//        { 
//        store_ = other.store_;
//        }
//    else               
//        { 
//        Error("Not yet implemented");
//        store_ = make_shared<ITDense<Real>>(is_);
//        applyFunc<Reshape>(store_,other.store_,{P,other.is_});
//        }
//    }


//class IsScalar
//    {
//    bool value_ = false;
//    public:
//
//    explicit operator bool() const { return value_; }
//
//    template <typename T>
//    NewData
//    operator()(const ITScalar<T>& d)
//        {
//        value_ = true;
//        return NewData();
//        }
//
//    template <class T>
//    NewData
//    operator()(const T& d)
//        {
//        value_ = false;
//        return NewData();
//        }
//    };

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

class Contract
    {
    const Label &Lind_,
                &Rind_;

    const IndexSet &Lis_,
                   &Ris_;

    //New IndexSet
    IndexSet Nis_;
    Real scalefac_ = -1;

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

    ITResult
    operator()(const ITDense<Real>& a1,
               const ITDense<Real>& a2);

    ITResult
    operator()(const ITDense<Real>& d,
               const ITCombiner& C)
        {
        auto res = combine(d,Lis_,Ris_);
        if(!res) return ITResult::None;
        else     return std::move(res);
        }
    ITResult
    operator()(const ITCombiner& C,
               const ITDense<Real>& d)
        { 
        auto res = combine(d,Ris_,Lis_);
        if(!res) return ITResult::AssignPointer;
        else     return std::move(res);
        }

    ITResult
    operator()(const ITDiag<Real>& d,
               const ITDense<Real>& t)
        {
        return diagDense(d,Lis_,Lind_,t,Ris_,Rind_);
        }
    ITResult
    operator()(const ITDense<Real>& t,
               const ITDiag<Real>& d)
        { 
        return diagDense(d,Ris_,Rind_,t,Lis_,Lind_);
        }

    //ITResult
    //operator()(const ITDense<Real>& a1,
    //           const ITDense<Complex>& a2) const
    //    {
    //    ITDense<Complex> c1(a1);
    //    return operator()(c1,a2);
    //    }

    //ITResult
    //operator()(const ITDense<Complex>& a1,
    //           const ITDense<Real>& a2) const
    //    {
    //    ITDense<Complex> c2(a2);
    //    return operator()(a1,c2);
    //    }


    //template <typename T1, typename T2>
    //ITResult
    //operator()(const ITDense<T1>& a1,
    //           const ITDense<T2>& a2) const
    //    {
    //    using product_type = decltype(::std::declval<T1>() * ::std::declval<T2>());
    //    //static const auto One = product_type(1.),
    //    //                  Zero = product_type(0.);
    //    auto res = new ITDense<product_type>();
    //    //TODO:
    //    Error("Contract not implemented for tensors of different element types.");
    //    //btas::contract(One,a1.t_,Lind_,a2.t_,Rind_,Zero,res->t_,Nind_);
    //    return ITResult(res);
    //    }

    //template <typename T1, typename T2>
    //ITResult
    //operator()(const T1& a1,const T2& a2) const
    //    {
    //    Error("Contract not implemented for this case");
    //    return ITResult();
    //    }

    private:

    NewData
    combine(const ITDense<Real>& d,
            const IndexSet& dis,
            const IndexSet& Cis);

    ITResult
    diagDense(const ITDiag<Real>& d,
              const IndexSet& dis,
              const Label& dind,
              const ITDense<Real>& t,
              const IndexSet& tis,
              const Label& tind);
 
    };

ITResult Contract::
operator()(const ITDense<Real>& a1,
           const ITDense<Real>& a2)
    {
    contractIS(Lis_,Lind_,Ris_,Rind_,Nis_,true);
    
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

    //PRI(Lind);
    //PRI(Rind);
    //PRI(Nind);

    auto res = make_newdata<ITDense<Real>>(area(Nis_),0.);
    auto t1 = make_tensorref(a1.data.data(),Lis_),
         t2 = make_tensorref(a2.data.data(),Ris_);
    auto tr = make_tensorref(res->data.data(),Nis_);
    contractloop(t1,Lind_,t2,Rind_,tr,Nind);
    scalefac_ = 0;
    for(auto elt : res->data)
        {
        scalefac_ += elt*elt;
        }
    scalefac_ = std::sqrt(scalefac_);
    if(scalefac_ != 0)
        {
        for(auto& elt : res->data)
            {
            elt /= scalefac_;
            }
        }
    return move(res);
    }

ITResult Contract::
diagDense(const ITDiag<Real>& d,
          const IndexSet& dis,
          const Label& dind,
          const ITDense<Real>& t,
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
        auto res = make_newdata<ITDense<Real>>(area(Nis_),0.);
        auto *pr = res->data.data();
        const auto *pt = t.data.data();

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
            auto* pd = d.data.data();
            assert(d.data.size() == dsize);
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
        return move(res);
        }
    else
        {
        //all of t's indices contracted with d
        //result will be diagonal
        if(d_ustride == 0) //all of d's inds contracted
            {
            // o scalar if all of d's inds contracted also
            Real val = 0;
            const auto *pt = t.data.data();
            if(d.allSame())
                {
                for(size_t J = 0; J < dsize; ++J)
                    val += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(dsize == d.data.size());
                auto *pd = d.data.data();
                for(size_t J = 0; J < dsize; ++J)
                    val += pd[J]*pt[J*t_cstride];
                }
            auto res = make_newdata<ITDiag<Real>>(val);
            return move(res);
            }
        else //some of d's inds uncontracted
            {
            // o element-wise product of d's data and t's diagonal
            auto res = make_newdata<ITDiag<Real>>(dsize,0.);
            auto *pr = res->data.data();
            const auto *pt = t.data.data();
            if(d.allSame())
                {
                for(size_t J = 0; J < dsize; ++J)
                    pr[J] += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(dsize == d.data.size());
                auto *pd = d.data.data();
                for(size_t J = 0; J < dsize; ++J)
                    pr[J] += pd[J]*pt[J*t_cstride];
                }
            return move(res);
            }
        }
    Error("Case not handled");
    return ITResult();
    }

NewData Contract::
combine(const ITDense<Real>& d,
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
        return NewData();
        }
    else
        {
        //dis doesn't have cind, replace
        //Cis[1], Cis[2], ... with cind
        //may need to reshape
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
            newind.reserve(dis.r()-Cis.r()+1);
            for(int j = 0; j < J1; ++j) 
                newind.push_back(dis[j]);
            newind.push_back(cind);
            for(int j = J1+Cis.r()-1; j < dis.r(); ++j) 
                newind.push_back(dis[j]);
            Nis_ = IndexSet(move(newind));
            return NewData();
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
            auto res = make_newdata<ITDense<Real>>(area(Nis_));
            auto td = make_tensorref(d.data.data(),dis);
            auto tr = make_tensorref(res->data.data(),rr);
            reshape(td,P,tr);
            return move(res);
            }
        }
    return NewData();
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

    auto C = applyFunc<Contract>(store_,other.store_,{Lis,Lind,Ris,Rind});

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

ITensor& ITensor::
operator/=(Real fac) 
    { 
    scale_ /= fac; 
    return *this; 
    }

ITensor& ITensor::
operator*=(Complex z)
    {
    if(z.imag() == 0) return operator*=(z.real());
    solo();
    applyFunc<MultComplex>(store_,{z});
    return *this;
    }

ITensor& ITensor::
operator/=(Complex z)
    {
    return operator*=(1./z);
    }

ITensor ITensor::
operator-() const 
    { 
    ITensor T(*this); 
    T.scale_ *= -1; 
    return T; 
    }

ITensor& ITensor::
noprime(IndexType type) 
    { 
    is_.noprime(type); 
    return *this; 
    }

//Set primeLevel of Index I to zero
ITensor& ITensor::
noprime(const Index& I) 
    { 
    is_.noprime(I); 
    return *this; 
    }

//Increase primeLevel of Indices by 1 (or optional amount inc)
ITensor& ITensor::
prime(int inc) 
    { 
    prime(All,inc); 
    return *this;
    }

//Increase primeLevel of Indices by 1 (or optional amount inc)
ITensor& ITensor::
prime(IndexType type, int inc) 
    { 
    is_.prime(type,inc); 
    return *this; 
    }

//Increase primeLevel of Index I by 1 (or optional amount inc)
ITensor& ITensor::
prime(const Index& I, int inc) 
    { 
    is_.prime(I,inc); 
    return *this; 
    }

//Change all Indices having primeLevel plevold to have primeLevel plevnew
ITensor& ITensor::
mapprime(int plevold, int plevnew, IndexType type)
    { 
    is_.mapprime(plevold,plevnew,type); 
    return *this; 
    }

ITensor& ITensor::
operator+=(const ITensor& other)
    {
    if(!*this) Error("Calling += on default constructed ITensor");
    if(!other) Error("Right-hand-side of += is default constructed");
    if(this == &other) return operator*=(2.);
    if(this->scale_.isZero()) return operator=(other);

    PlusEQ::permutation P(is_.size());
    try {
        detail::calc_permutation(is_,other.is_,P);
        }
    catch(const ITError& e)
        {
        Print(*this);
        Print(other);
        Error("ITensor::operator+=: different Index structure");
        }

    Real scalefac = 1;
    if(scale_.magnitudeLessThan(other.scale_)) 
        {
        this->scaleTo(other.scale_); 
        }
    else
        {
        scalefac = (other.scale_/scale_).real();
        }

    solo();
    
    if(isTrivial(P))
        {
        applyFunc<PlusEQ>(store_,other.store_,{scalefac});
        }
    else
        {
        applyFunc<PlusEQ>(store_,other.store_,{P,is_,other.is_,scalefac});
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

ITensor& ITensor::
fill(Complex z)
    {
    if(!(*this)) return *this;
    solo();
    scale_ = LogNumber(1.);
    if(z.imag() == 0)
        applyFunc<FillReal>(store_,{z.real()});
    else
        applyFunc<FillCplx>(store_,{z});
    return *this;
    }

class MultReal
    {
    Real r_;
    public:
    MultReal(Real r)
        : r_(r)
        { }

    template<typename T>
    ITResult
    operator()(ITDense<T>& d) const
        {
        //TODO: use BLAS algorithm?
        for(auto& elt : d.data)
            elt *= r_;
        return ITResult();
        }

    template<typename T>
    ITResult
    operator()(ITDiag<T>& d) const
        {
        d.val *= r_;
        //TODO: use BLAS algorithm
        for(auto& elt : d.data)
            elt *= r_;
        return ITResult();
        }

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("MultReal not implemented for ITData type."); return ITResult(); }
    };

void ITensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    solo();
    scale_ /= newscale;
    applyFunc<MultReal>(store_,{scale_.real0()});
    scale_ = newscale;
    }


void ITensor::
solo()
	{
    if(!store_.unique()) store_ = store_->clone();
    }

class NormNoScale
    {
    Real nrm_;
    public:

    NormNoScale() : nrm_(0) { }

    operator Real() const { return nrm_; }

    template<typename T>
    ITResult
    operator()(const ITDense<T>& d) { return calc(d); }
    template<typename T>
    ITResult
    operator()(const ITDiag<T>& d) { return calc(d); }

    template<typename T>
    ITResult
    calc(const T& d)
        {
        for(const auto& elt : d.data)
            {
            nrm_ += std::norm(elt);
            }
        nrm_ = std::sqrt(nrm_);
        return ITResult();
        }
    };

void ITensor::
scaleOutNorm()
    {
    Real f = applyFunc<NormNoScale>(store_);
    //If norm already 1 return so
    //we don't have to call solo()
    if(fabs(f-1) < 1E-12) return;
    if(f == 0)
        {
        scale_ = LogNumber(1.);
        return;
        }

    solo();
    applyFunc<MultReal>(store_,{1./f});
    scale_ *= f;
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
        if(t) applyFunc<PrintIT>(t.data(),{s,t.scale(),t.inds()});
        else           s << " (default constructed)}\n";
        }
    else
        {
        s << "}";
        }
    return s;
    }



ITensor
randomize(ITensor T, const Args& args)
    {
    T.generate(detail::quickran);
    return T;
    }


Real
norm(const ITensor& T)
    {
#ifdef DEBUG
    if(!T) Error("ITensor is default initialized");
#endif
    return T.scale().real0() *
           applyFunc<NormNoScale>(T.data());
    }

//Possible optimization:
// this will get called a lot,
// optimize by requiring all ITData subtypes to
// have a virtual method bool isComplex() ?
ITensor
conj(const ITensor& T)
    {
    if(isComplex(T))
        {
        auto Tc = T;
        Tc.apply([](auto z) { return std::conj(z); });
        return Tc;
        }
    return T;
    }


class CheckComplex
    {
    bool isComplex_;
    public:

    CheckComplex() : isComplex_(false) { }

    operator bool() const { return isComplex_; }

    NewData
    operator()(const ITDense<Real>& d) { isComplex_ = false; return NewData(); }
    NewData
    operator()(const ITDense<Complex>& d) { isComplex_ = true; return NewData(); }

    template <class T>
    NewData
    operator()(const T& d)
        {
        Error("CheckComplex not implemented for data type.");
        return NewData();
        }
    };

bool
isComplex(const ITensor& t)
    {
    return applyFunc<CheckComplex>(t.data());
    }

class SumEls
    {
    Complex sum_;
    public:

    SumEls() : sum_(0) { }

    operator Complex() const { return sum_; }

    template <class T>
    NewData
    operator()(const ITDense<T>& d) 
        { 
        for(const auto& elt : d.data)
            sum_ += elt;
        return NewData();
        }

    template <class T>
    NewData
    operator()(const ITDiag<T>& d) 
        { 
        for(const auto& elt : d.data)
            sum_ += elt;
        return NewData();
        }
    };
Real
sumels(const ITensor& t)
    {
    auto z = Complex(applyFunc<SumEls>(t.data()));
    if(z.imag() != 0) Error("ITensor has non-zero imaginary part");
    return t.scale().real0()*z.real();
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
    return ITensor(IndexSet(std::move(inds)),make_newdata<ITCombiner>(),{1.0});
    }

ITensor
delta(const Index& i1, const Index& i2)
    {
#ifdef DEBUG
    if(i1.m() != i2.m()) Error("delta: indices must have same dimension");
#endif
    return ITensor({i1,i2},make_newdata<ITCombiner>(),{1.0});
    }

};
