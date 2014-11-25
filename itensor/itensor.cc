//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor.h"
#include "contract.h"

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
ITensor()  
    { }


ITensor::
ITensor(const Index& i1) 
    :
    is_(i1),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(i1.m()))
	{ 
    }


ITensor::
ITensor(const Index& i1,const Index& i2) 
    :
    is_(i1,i2),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(i1.m(),i2.m()))
	{ 
    }
    

ITensor::
ITensor(Real val) 
    :
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(1))
    { 
    fill(val);
    }

ITensor::
ITensor(Complex val) 
    :
    scale_(1.),
    store_(std::make_shared<ITDense<Complex>>(1))
    { 
    fill(val);
    }

ITensor::
ITensor(IndexSet<Index>&& iset,
        NewData nd,
        LogNumber scale)
    :
    is_(iset),
    scale_(scale),
    store_(std::move(nd))
    {
    }

ITensor::
ITensor(const Index& i1,
        const VectorRef& V)
    :
    is_(i1),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(is_,V.begin(),V.end()))
	{ 
    }

ITensor::
ITensor(const Index& i1,
        const Index& i2,
        const VectorRef& V)
    :
    is_(i1,i2),
    scale_(1.),
    store_(std::make_shared<ITDiag<Real>>(V.begin(),V.end()))
	{ 
#ifdef DEBUG
    if(V.Length() != std::min(i1.m(),i2.m()))
        Error("Wrong size of data in diagonal ITensor constructor");
#endif
    }

ITensor::
ITensor(const IndexSet<Index>& is)
    :
    is_(is),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(is_))
	{ }

ITensor::
ITensor(const IndexSet<Index>& is,
        const VectorRef& v)
    :
    is_(is),
    scale_(1.),
    store_(std::make_shared<ITDense<Real>>(is_,v.begin(),v.end()))
	{ }

struct CopyElems
    {
    CopyElems() { }

    template<typename T>
    NewData
    operator()(T& d1,
               const T& d2)
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
ITensor(const IndexSet<Index>& is,
        const ITensor& t)
    :
    is_(is),
    scale_(t.scale_),
    store_(std::make_shared<ITDense<Real>>(is_))
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
        store_ = make_newdata<ITDense<Real>>(is_,vref.begin(),vref.end());
        }
    else
        {
        VectorRef vref = M.TreatAsVector(); 
        store_ = make_newdata<ITDense<Real>>(is_,vref.begin(),vref.end());
        }
    }


//class Reshape
//    {
//    const Permutation& P_;
//    const IndexSet<Index>& is_;
//    public:
//    Reshape(const Permutation& P,
//            const IndexSet<Index>& is)
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
//ITensor(const IndexSet<Index>& is,
//        const ITensor& t,
//        const Permutation& P)
//    :
//    is_(is_),
//    scale_(t.scale_)
//    {
//    if(P.isTrivial()) 
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


ITensor& ITensor::
operator*=(const ITensor& other)
    {
    if(!(*this) || !other)
        Error("Default constructed ITensor in product");

    if(this == &other)
        {
        return operator=( ITensor(sqr(this->norm())) );
        }

    //Check if other is a scalar
    if(other.inds().r() == 0)
        {
        return operator*=(other.cplx());
        }
    //Check if this is a scalar
    if(this->r() == 0)
        {
        return operator=(other*cplx());
        }

    std::vector<bool> contL(is_.r(),false),
                      contR(other.is_.r(),false);

    //Set Lind, Rind to zero. Special value 0 marks
    //uncontracted indices. Later will assign unique numbers
    //to these entries in Lind and Rind
    Labels Lind(is_.rn(),0),
           Rind(other.is_.rn(),0);

    //Count number of contracted indices,
    //set corresponding entries of Lind, Rind
    //to 1,2,...,ncont
    int ncont = 0;
    for(int i = 0; i < is_.rn(); ++i)
    for(int j = 0; j < other.is_.rn(); ++j)
        if(is_[i] == other.is_[j])
            {
            contL[i] = true;
            contR[j] = true;

            ++ncont;
            Lind[i] = ncont;
            Rind[j] = ncont;

            break;
            }

    //Finish making contL, contR for m==1 indices
    int ncont_all = ncont;
    for(int i = is_.rn(); i < is_.r(); ++i)
    for(int j = other.is_.rn(); j < other.is_.r(); ++j)
        {
        if(is_[i] == other.is_[j])
            {
            ++ncont_all;
            contL[i] = true;
            contR[j] = true;
            break;
            }
        }

    //nuniq is total number of unique, uncontracted indices
    //(nuniq all includes m==1 indices)
    int nuniq = is_.rn()+other.is_.rn()-2*ncont;
    int nuniq_all = is_.r()+other.is_.r()-2*ncont_all;

    //container in which we will accumulate the new indices
    IndexSet<Index>::storage newind;
    newind.reserve(nuniq_all);

    //Go through and assign uncontracted entries of Lind,Rind
    //the integers ncont+1,ncont+2,...
    //Simultaneously fill newind (keeping count "ni")
    int uu = ncont;
    for(int j = 0; j < is_.rn(); ++j)
        {
        if(!contL[j]) 
            {
            Lind[j] = ++uu;
            newind.push_back(is_[j]);
            }
        }
    for(int j = 0; j < other.is_.rn(); ++j)
        {
        if(!contR[j]) 
            {
            Rind[j] = ++uu;
            newind.push_back(other.is_[j]);
            }
        }

    //Finish filling up newind with m==1 indices
    for(int j = is_.rn(); j < is_.r(); ++j)
        {
        if(!contL[j]) newind.push_back(is_[j]);
        }
    for(int j = other.is_.rn(); j < other.is_.r(); ++j)
        {
        if(!contR[j]) newind.push_back(other.is_[j]);
        }

    IndexSet<Index> new_index(std::move(newind));
#ifdef DEBUG
    if(new_index.rn() != nuniq) Error("new_index size not equal to nuniq");
#endif

    Labels Pind(nuniq);
    for(int i = 0; i < new_index.r(); ++i)
        {
        int j = hasindex(is_,new_index[i]);
        if(j)
            {
            Pind[i] = Lind[j];
            }
        else
            {
            j = hasindex(other.is_,new_index[i]);
            Pind[i] = Rind[j];
            }
        }

    //println(this->is_);
    //cout << "Lind = {";
    //for(auto x : Lind)
    //    {
    //    cout << x << ",";
    //    }
    //cout << "}" << endl;
    //println(other.is_);
    //cout << "Rind = {";
    //for(auto x : Rind)
    //    {
    //    cout << x << ",";
    //    }
    //cout << "}" << endl;
    //println(new_index);
    //cout << "Pind = {";
    //for(auto x : Pind)
    //    {
    //    cout << x << ",";
    //    }
    //cout << "}" << endl;
    //exit(0);

    applyFunc<Contract>(store_,other.store_,{Lind,Rind,Pind,new_index});

    is_.swap(new_index);

    scale_ *= other.scale_;

    scaleOutNorm();

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

    if(is_ != other.is_)
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

    applyFunc<PlusEQ>(store_,other.store_,{scalefac});

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

    if(f != 0)
        {
        solo();
        applyFunc<MultReal>(store_,{1./f});
        scale_ *= f;
        }
    else //norm == zero
        {
        scale_ = LogNumber(1.);
        }
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

Real ITensor::
norm() const 
    {
#ifdef DEBUG
    if(!*this) Error("ITensor is default initialized");
#endif
    return scale_.real0() *
           applyFunc<NormNoScale>(store_);
    }


ostream& 
operator<<(ostream & s, const ITensor& t)
    {
    s << "ITensor r = " << t.r() << ": ";
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
    return s;
    }

ITensor
randIT(ITensor T, const Args& args)
    {
    T.generate(Global::random);
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
    operator()(const T& d) 
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

};
