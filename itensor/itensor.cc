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
    auto ncont = computeAnnotations(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind);
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
fill(Real r)
    {
    if(!(*this)) return *this;
    solo();
    scale_ = LogNumber(1.);
    applyFunc<FillReal>(store_,{r});
    return *this;
    }

ITensor& ITensor::
fill(Complex z)
    {
    if(!(*this)) return *this;
    solo();
    scale_ = LogNumber(1.);
    applyFunc<FillCplx>(store_,{z});
    return *this;
    }

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


Real
quickran()
    {
    static auto seed = (std::time(NULL) + getpid());
    int im = 134456;
    int ia = 8121;
    int ic = 28411;
    Real scale = 1.0 / im;
    seed = (seed*ia+ic)%im;
    return Real(seed) * scale;
    }

ITensor
randIT(ITensor T, const Args& args)
    {
    T.generate(quickran);
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
