#include "itensor/itensor_interface.h"

namespace itensor {

Cplx ITensor::
cplx() const
    {
    if(inds().r() != 0)
        {
        Error(format("Wrong number of IndexVals passed to real/cplx (expected %d, got 0)",inds().r()));
        }
    constexpr size_t size = 0;
    auto inds = IntArray(size);
    auto z = itensor::doTask(GetElt{is_,inds},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale_.real0(); 
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in cplx(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in cplx(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

void ITensor::
set(Cplx val)
    {
    if(0 != size_t(inds().r())) 
        {
        Error(format("Wrong number of IndexVals passed to set (expected %d, got 0)",
                     inds().r()));
        }
    auto inds = IntArray(0,1);
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.)
        {
        doTask(SetElt<Real>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{val,is_,inds},store_);
        }
    }

void ITensor::
set(std::vector<int> const& ints,
    Cplx val)
    {
    auto size = ints.size();
    if(size != size_t(inds().r())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ints) println(iv);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().r(),size));
        }
    auto inds = IntArray(is_.r(),0);
    for(auto i : range(size))
        inds[i] = ints[i]-1;
    //TODO: if !store_ and !is_real, call allocCplx instead
    //detail::permute_map(is_,ivals,inds,
    //                    [](IndexVal const& iv) { return iv.val-1; });
    if(!store_) detail::allocReal(*this,inds); 
    scaleTo(1.);
    if(val.imag()==0.0)
        {
        doTask(SetElt<Real>{val.real(),is_,inds},store_);
        }
    else
        {
        doTask(SetElt<Cplx>{val,is_,inds},store_);
        }
    }

ITensor& ITensor::
conj()
    {
    doTask(Conj{},store_);
    return *this;
    }

ITensor& ITensor::
dag()
    {
    Error("dag not implemented");
    return *this;
    }

ITensor& ITensor::
takeReal()
    {
    doTask(TakeReal{},store_);
    return *this;
    }

ITensor& ITensor::
takeImag()
    {
    doTask(TakeImag{},store_);
    return *this;
    }


ITensor& ITensor::
fill(Cplx z)
    {
    if(!store_) 
        {
        if(is_) detail::allocReal(*this);
        else Error("Can't fill default-constructed tensor");
        }
    IF_USESCALE(scale_ = scale_type(1.);)
    if(z.imag() == 0)
        doTask(Fill<Real>{z.real()},store_);
    else
        doTask(Fill<Cplx>{z},store_);
    return *this;
    }

#ifdef USESCALE
void ITensor::
scaleTo(scale_type const& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    doTask(Mult<Real>{scale_.real0()},store_);
    scale_ = newscale;
    }

void inline ITensor::
scaleTo(Real newscale) { scaleTo(LogNum{newscale}); }
#endif

void ITensor::
swap(ITensor & other)
    {
    is_.swap(other.is_);
    store_.swap(other.store_);
    IF_USESCALE(scale_.swap(other.scale_);)
    }

Index
commonIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return Index{};
    }

Index
uniqueIndex(ITensor const& A, 
            ITensor const& B, 
            IndexType t)
    {
    for(auto& I : A.inds())
        if( (t == All || I.type() == t)
         && !hasindex(B.inds(),I) ) 
            {
            return I;
            }
    return Index{};
    }

ITensor
swapPrime(ITensor T, 
          int plev1, 
          int plev2,
          IndexType type)
    { 
    int tempLevel = 99999;
#ifdef DEBUG
    for(auto& I : T.inds())
        {
        if(I.primeLevel() == tempLevel) 
            {
            println("tempLevel = ",tempLevel);
            Error("swapPrime fails if an index has primeLevel==tempLevel");
            }
        }
#endif
    T.mapprime(plev1,tempLevel,type);
    T.mapprime(plev2,plev1,type);
    T.mapprime(tempLevel,plev2,type);
    return T; 
    }

Real
norm(ITensor const& T)
    {
#ifdef DEBUG
    if(!T) Error("Default initialized tensor in norm(ITensorT)");
#endif

#ifndef USESCALE
    return doTask(NormNoScale{},T.store());
#else
    auto fac = std::fabs(T.scale().real0());
    return fac * doTask(NormNoScale{},T.store());
#endif
    }

void
randomize(ITensor & T, Args const& args)
    {
    if(!T.store()) detail::allocReal(T);
#ifdef DEBUG
    if(!T) Error("default initialized tensor in randomize");
#endif
    auto cplx = args.getBool("Complex",false);
    if(cplx) T.generate(detail::quickranCplx);
    else     T.generate(detail::quickran);
    }

void ITensor::
write(std::ostream& s) const
    {
    itensor::write(s,inds());
    itensor::write(s,scale());
    auto type = StorageType::Null;
    if(store()) 
        {
        type = doTask(StorageType{},store());
        }
    itensor::write(s,type);
    if(store()) 
        {
        doTask(Write{s},store());
        }
    }

ITensor
multSiteOps(ITensor A, ITensor const& B) 
    {
    A.prime(Site);
    A *= B;
    A.mapprime(2,1,Site);
    return A;
    }

bool
hasindex(ITensor const& T, Index const& I)
    {
    return detail::contains(T.inds(),I);
    }

Index
findtype(ITensor const& T, IndexType type)
    {
    for(auto& i : T.inds())
        if(i.type()==type) return i;
    return Index{};
    }

bool
isComplex(ITensor const& T)
    {
    return doTask(CheckComplex{},T.store());
    }

bool
isReal(ITensor const& T)
    {
    return not isComplex(T);
    }

void ITensor::
read(std::istream& s)
    {
    itensor::read(s,is_);
    LogNum scale;
    itensor::read(s,scale);
    IF_USESCALE(scale_ = scale;)
    auto type = StorageType::Null;
    itensor::read(s,type);
    if(type==StorageType::Null) { /*intentionally left blank*/  }
    else if(type==StorageType::DenseReal) { store_ = readType<DenseReal>(s); }
    else if(type==StorageType::DenseCplx) { store_ = readType<DenseCplx>(s); }
    else if(type==StorageType::Combiner) { store_ = readType<Combiner>(s); }
    else if(type==StorageType::DiagReal) { store_ = readType<Diag<Real>>(s); }
    else if(type==StorageType::DiagCplx) { store_ = readType<Diag<Cplx>>(s); }
    //else if(type==StorageType::QDenseReal) { store_ = readType<QDense<Real>>(s); }
    //else if(type==StorageType::QDenseCplx) { store_ = readType<QDense<Cplx>>(s); }
    //else if(type==StorageType::QDiagReal) { store_ = readType<QDiag<Real>>(s); }
    //else if(type==StorageType::QDiagCplx) { store_ = readType<QDiag<Cplx>>(s); }
    //else if(type==StorageType::QCombiner) { store_ = readType<QCombiner>(s); }
    //else if(type==StorageType::ScalarReal) { store_ = readType<ScalarReal>(s); }
    //else if(type==StorageType::ScalarCplx) { store_ = readType<ScalarCplx>(s); }
    else
        {
        Error("Unrecognized type when reading tensor from istream");
        }
    }

} //namespace itensor
