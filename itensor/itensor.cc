//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/print_macro.h"
#include "itensor/util/range.h"
#include "itensor/util/safe_ptr.h"
#include "itensor/itensor.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/contract.h"

using std::array;
using std::ostream;
using std::vector;

namespace itensor {

//
// ITensor Constructors
//

ITensor::
ITensor(std::vector<index_type> const& inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(std::initializer_list<index_type> inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(indexset_type const& is)
  : is_(is)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(indexset_type iset,
        storage_ptr&& pdat,
        LogNum const& scale)
    :
    is_(std::move(iset)),
    store_(std::move(pdat))
    { 
    IF_USESCALE(scale_ = scale;)
    }

    
ITensor::
ITensor(Cplx val) 
    { 
    //IF_USESCALE(scale_ = LogNum(1.);)
    //if(val.imag() == 0)
    //    {
    //    store_ = newITData<ScalarReal>(val.real());
    //    }
    //else
    //    {
    //    store_ = newITData<ScalarCplx>(val);
    //    }
    ////if(val.imag() == 0)
    ////    store_ = newITData<Diag<Real>>(1,val.real());
    ////else
    ////    store_ = newITData<Diag<Cplx>>(1,val);
    }


//template<>
//void ITensor::
//scaleOutNorm()
//    {
//    auto nrm = doTask(NormNoScale<Index>{is_},store_);
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
set(std::vector<IndexVal> const& ivals,
    Cplx val)
    {
    auto size = ivals.size();
    if(size != size_t(inds().r())) 
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ivals) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to set (expected %d, got %d)",
                     inds().r(),size));
        }
    auto inds = IntArray(is_.r(),0);
    detail::permute_map(is_,ivals,inds,
                        [](IndexVal const& iv) { return iv.val-1; });
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
dag() { return conj(); }

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


ITensor 
random(ITensor T, const Args& args)
    {
    randomize(T,args);
    return T;
    }

Real 
sumels(ITensor const& t)
    {
    auto z = sumelsC(t);
    if(z.imag() != 0) Error("ITensor has non-zero imaginary part, use sumelsC");
    return z.real();
    }

Cplx 
sumelsC(ITensor const& t)
    {
    auto z = doTask(SumEls{t.inds()},t.store());
#ifndef USESCALE
    return z;
#else
    return t.scale().real0()*z;
#endif
    }

ostream& 
operator<<(ostream & s, ITensor const& t)
    {
    s << "ITensor r=" << t.r() << ": " << t.inds() << "\n";
    if(!t.store()) 
        {
        s << "{Zero / Not yet allocated}\n";
        }
    else
        {
        //Checking whether std::ios::floatfield is set enables 
        //printing the contents of an ITensor when using the printf
        //format string %f (or another float-related format string)
        bool ff_set = (std::ios::floatfield & s.flags()) != 0;
        bool print_data = (ff_set || Global::printdat());
        doTask(PrintIT{s,t.scale(),t.inds(),print_data},t.store());
        }
    return s;
    }

ITensor
matrixTensor(Matrix&& M, Index const& i1, Index const& i2)
    {
    auto res = ITensor({i1,i2},DenseReal{std::move(M.storage())});
    M.clear();
    return res;
    }

ITensor
matrixTensor(Matrix const& M, Index const& i1, Index const& i2)
    {
    return matrixTensor(Matrix(M),i1,i2);
    }

ITensor
matrixTensor(CMatrix&& M, Index const& i1, Index const& i2)
    {
    bool isReal = true;
    for(auto& el : M)
    if(std::fabs(el.imag()) > 1E-14)
        {
        isReal = false;
        break;
        }
    ITensor res;
    if(isReal)
        {
        auto store = vector<Real>(M.size());
        for(auto n : range(M.size())) store[n] = M.store()[n].real();
        res = ITensor({i1,i2},DenseReal{std::move(store)});
        }
    else
        {
        res = ITensor({i1,i2},DenseCplx{std::move(M.storage())});
        }
    M.clear();
    return res;
    }

ITensor
matrixTensor(CMatrix const& M, Index const& i1, Index const& i2)
    {
    return matrixTensor(CMatrix(M),i1,i2);
    }


ITensor
combiner(std::vector<Index> inds, Args const& args)
    {
    if(inds.empty()) Error("No indices passed to combiner");
    long rm = 1;
    for(const auto& i : inds)rm *= i.m();
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
    return ITensor(IndexSet(std::move(inds)),Combiner{});
    }

struct IsCombiner
    {
    template<typename D>
    bool 
    operator()(D const& d) { return false; }
    bool
    operator()(Combiner const& d) { return true; }
    };

Index
combinedIndex(ITensor const& C)
    {
#ifdef DEBUG
    auto iscombiner = applyFunc(IsCombiner{},C.store());
    if(not iscombiner)
        {
        throw ITError("Called combinedIndex on ITensor that is not a combiner");
        }
#endif
    return C.inds().front();
    }

ITensor
randomTensor(IndexSet const& inds)
    {
    return random(ITensor{inds});
    }

QN
div(ITensor const& T)
    {
    Error("div not yet implemented");
    return QN();
    }

} //namespace itensor
