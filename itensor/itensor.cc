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
using std::move;

namespace itensor {

//
// ITensor Constructors
//

ITensor::
ITensor(std::vector<Index> const& inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(std::initializer_list<Index> inds)
  : is_(inds)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(IndexSet const& is)
  : is_(is)
    { 
    IF_USESCALE(scale_ = LogNum(1.);)
    }

ITensor::
ITensor(IndexSet iset,
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
    IF_USESCALE(scale_ = LogNum(1.);)
    if(val.imag() == 0)
        {
        store_ = newITData<ScalarReal>(val.real());
        }
    else
        {
        store_ = newITData<ScalarCplx>(val);
        }
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
eltC() const
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
        println("too big for real in eltC(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in eltC(...)");
        return Cplx(0.,0.);
        }
    return Cplx(NAN,NAN);
#endif
    }

Cplx ITensor::
eltC(std::vector<IndexVal> const& ivs) const
    {
    if(!store()) Error("tensor storage unallocated");

    auto size = ivs.size();
    if(size != size_t(inds().r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",inds(),"\n");
        println("---------------------------------------------");
        println("Indices provided = ");
        for(auto& iv : ivs) println(iv.index);
        println("---------------------------------------------");
        Error(format("Wrong number of IndexVals passed to real/cplx (expected %d, got %d)",inds().r(),size));
        }

    auto ints = IntArray(size);
    detail::permute_map(inds(),ivs,ints,
                [](IndexVal const& iv) { return iv.val-1; });
    auto z = itensor::doTask(GetElt{inds(),ints},store_);
#ifndef USESCALE
    return z;
#else
    try {
        return z*scale().real0();
        }
    catch(TooBigForReal const& e)
        {
        println("too big for real in eltC(...), scale = ",scale());
        throw e;
        }
    catch(TooSmallForReal const&)
        {
        println("warning: too small for real in eltC(...)");
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
dag() 
    { 
    if(hasQNs(is_)) is_.dag();
    return conj(); 
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
randomize(Args const& args)
    {
    if(!this->store()) detail::allocReal(*this);
#ifdef DEBUG
    if(!(*this)) Error("default initialized tensor in randomize");
#endif
    auto cplx = args.getBool("Complex",false);
    if(cplx) this->generate(detail::quickranCplx);
    else     this->generate(detail::quickran);
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

void ITensor::
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
            TagSet const& ts)
    {
    for(auto& I : A.inds())
        if( (hasTags(ts,TagSet(All)) || hasTags(I,ts))
         && hasIndex(B.inds(),I) ) 
            {
            return I;
            }
    return Index();
    }

Index
uniqueIndex(ITensor const& A, 
            ITensor const& B, 
            TagSet const& ts)
    {
    for(auto& I : A.inds())
        if( (hasTags(ts,TagSet(All)) || hasTags(I,ts))
         && !hasIndex(B.inds(),I) ) 
            {
            return I;
            }
    return Index();
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
    Global::warnDeprecated("randomize(ITensor,args) is deprecated in favor of .randomize(args)");
    T.randomize(args);
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
    else if(type==StorageType::QDenseReal) { store_ = readType<QDense<Real>>(s); }
    else if(type==StorageType::QDenseCplx) { store_ = readType<QDense<Cplx>>(s); }
    else if(type==StorageType::QDiagReal) { store_ = readType<QDiag<Real>>(s); }
    else if(type==StorageType::QDiagCplx) { store_ = readType<QDiag<Cplx>>(s); }
    else if(type==StorageType::QCombiner) { store_ = readType<QCombiner>(s); }
    else if(type==StorageType::ScalarReal) { store_ = readType<ScalarReal>(s); }
    else if(type==StorageType::ScalarCplx) { store_ = readType<ScalarCplx>(s); }
    else
        {
        Error("Unrecognized type when reading tensor from istream");
        }
    }

namespace detail {

void
allocReal(ITensor& T)
    {
    if(hasQNs(T)) Error("Can't allocate quantum ITensor with undefined divergence");
    T.store() = newITData<DenseReal>(area(T.inds()),0);
    }

void
allocReal(ITensor& T, IntArray const& inds)
    {
    if(not hasQNs(T))
        {
        T.store() = newITData<DenseReal>(area(T.inds()),0);
        }
    else
        {
        QN div;
        for(auto i : range(T.inds()))
            {
            auto iv = (T.inds()[i])(1+inds[i]);
            div += iv.qn()*iv.index.dir();
            }
        T.store() = newITData<QDenseReal>(T.inds(),div);
        }
    }

void
allocCplx(ITensor& T)
    {
    if(hasQNs(T)) Error("Can't allocate quantum ITensor with undefined divergence");
    T.store() = newITData<DenseCplx>(area(T.inds()),0);
    }


void
checkArrows(IndexSet const& is1,
            IndexSet const& is2,
            bool shouldMatch = false)
    {
    if(hasQNs(is1) && hasQNs(is2))
        {
        for(auto I1 : is1)
        for(auto I2 : is2)
            {
            if(I1 == I2)
                {
                auto cond = shouldMatch ^ (I1.dir() == I2.dir());
                if(cond)
                    {
                    println("----------------------------------------");
                    println("IQIndexSet 1 = \n",is1);
                    println("----------------------------------------");
                    println("IQIndexSet 2 = \n",is2);
                    println("----------------------------------------");
                    printfln("Mismatched IQIndex from set 1 %s",I1);
                    printfln("Mismatched IQIndex from set 2 %s",I2);
                    Error("Mismatched IQIndex arrows");
                    }
                }
            }
        }
    }

void
checkSameDiv(ITensor const& T1,
             ITensor const& T2)
    {
    if(hasQNs(T1) && hasQNs(T2))
        {
        if(div(T1) != div(T2)) 
            {
            Error(format("div(T1)=%s must equal div(T2)=%s when adding T1+T2",div(T1),div(T2)));
            }
        }
    }

} //namespace detail

//TODO: implement proper Dense*QDense to avoid conversion cost
ITensor& ITensor::
operator*=(ITensor const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(L.r() == 0)
        {
        auto z = L.eltC();
        *this = R*z;
        return *this;
        }
    else if(R.r()==0)
        {
        auto z = R.eltC();
        *this *= z;
        return *this;
        }

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds());

    //TODO: create a proper doTask(Contract,Dense,QDense)
    auto hqL = hasQNs(L);
    auto hqR = hasQNs(R);
    auto Rdense = R;
    if(hqL && !hqR) L = removeQNs(L);
    else if(!hqL && hqR) Rdense = removeQNs(Rdense);

    auto C = doTask(Contract{L.inds(),Rdense.inds()},
                    L.store(),
                    Rdense.store());

#ifdef USESCALE
    L.scale_ *= Rdense.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;
#endif

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }

ITensor& ITensor::
permute(IndexSet const& iset)
    {
    auto& A = *this;
    auto Ais = A.inds();
    auto r = Ais.r();

    if(size_t(r) != size_t(iset.r()))
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",Ais,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n",iset,"\n");
        println("---------------------------------------------");
        Error(format("Wrong number of Indexes passed to permute (expected %d, got %d)",r,iset.r()));
        }

    // Get permutation
    auto P = Permutation(r);
    calcPerm(Ais,iset,P);
    if(isTrivial(P))
        {
        return A;
        }
    // If not trivial, use permutation to get new index set
    // This is necessary to preserve the proper arrow direction of IQIndex
    auto bind = RangeBuilderT<IndexSet>(r);
    for(auto i : range(r))
        {
        bind.setIndex(P.dest(i),Ais[i]);
        }
    auto Bis = bind.build();

    auto O = Order{P,Ais,Bis};
    if(A.store())
        {
        doTask(O, A.store());
        }

    A.is_.swap(Bis);

    return A;
    }

#ifndef USESCALE

////for Diag and QDiag
////QDense -> Dense
////QDiag  -> Dense
////Diag   -> Dense
//ITensor
//toDense(ITensor T)
//    {
//    if(not hasQNs(T)) return T;
//    if(T.store()) doTask(ToDense{T.inds()},T.store());
//    auto nis = T.inds();
//    nis.removeQNs();
//    return ITensor{move(nis),move(T.store()),T.scale()};
//    }

//TODO: make this use a RemoveQNs task type that does:
//QDense -> Dense
//QDiag  -> Diag
ITensor
removeQNs(ITensor T)
    {
    if(not hasQNs(T)) return T;
    if(T.store()) doTask(RemoveQNs{T.inds()},T.store());
    auto nis = T.inds();
    nis.removeQNs();
    return ITensor{move(nis),move(T.store()),T.scale()};
    }

ITensor& ITensor::
operator*=(Real r)
    {
    doTask(Mult<Real>{r},store_);
    return *this;
    }

ITensor& ITensor::
operator/=(Real r)
    {
    auto fac = 1./r;
    doTask(Mult<Real>{fac},store_);
    return *this;
    }

#endif

ITensor& ITensor::
operator*=(Cplx z)
    {
    if(z.imag() == 0) return operator*=(z.real());
    doTask(Mult<Cplx>{z},store_);
    return *this;
    }

//Non-contracting product
ITensor& ITensor::
operator/=(ITensor const& R)
    {
    auto& L = *this;

    if(!L || !R) Error("Default constructed ITensor in product");

    if(Global::checkArrows()) detail::checkArrows(L.inds(),R.inds(),true);

    auto C = doTask(NCProd{L.inds(),R.inds()},
                    L.store(),
                    R.store());

#ifdef USESCALE
    L.scale_ *= R.scale();
    if(!std::isnan(C.scalefac)) L.scale_ *= C.scalefac;
#endif

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    L.is_.swap(C.Nis);

    return L;
    }

void
daxpy(ITensor & L,
      ITensor const& R,
      Real alpha)
    {
    if(L.r() != R.r()) Error("ITensor::operator+=: different number of indices");

    using permutation = typename PlusEQ::permutation;

    auto P = permutation(L.inds().size());

    try {
        calcPerm(R.inds(),L.inds(),P);
        }
    catch(std::exception const& e)
        {
        println("L = ",L);
        println("R = ",R);
        Error("ITensoITensor::operator+=: different index structure");
        }

    if(Global::checkArrows()) 
        {
        auto shouldMatch = true;
        detail::checkArrows(L.inds(),R.inds(),shouldMatch);
        }

    if(!L.store()) Error("L not initialized in daxpy");

#ifdef DEBUG
    detail::checkSameDiv(L,R);
#endif

#ifdef USESCALE
    if(L.scale().magnitudeLessThan(R.scale())) 
        {
        L.scaleTo(R.scale()); 
        }
    else
        {
        alpha *= (R.scale()/L.scale()).real();
        }
#endif

    auto PEq = PlusEQ{P,L.inds(),R.inds(),alpha};
    doTask(PEq,L.store(),R.store());
    } 

ITensor& ITensor::
operator+=(ITensor const& R)
    {
    auto& L = *this;
    if(!L || !L.store()) { return (L=R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor += is default constructed");
    if(&L == &R) return operator*=(2.);

    daxpy(L,R,1.);
    
    return L;
    } 

ITensor& ITensor::
operator-=(ITensor const& R)
    {
    auto& L = *this;
    if(!L || !L.store()) { return (L = -R); } //special case when this (L) is not initialized
    if(!R) Error("Right-hand-side of ITensor -= is default constructed");
    if(&L == &R) 
        { 
        L *= 0.;
        return L;
        }
    daxpy(L,R,-1.);
    return L;
    }

ITensor
multSiteOps(ITensor A, ITensor const& B) 
    {
    A.prime("Site");
    A *= B;
    A.mapPrime(2,1,"Site");
    return A;
    }

bool
hasIndex(ITensor const& T, Index const& I)
    {
    return detail::contains(T.inds(),I);
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
    T.randomize(args);
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
    s << "ITensor ord=" << order(t) << ": "; 
    if(hasQNs(t)) 
        {
        if(t.r() > 0) s << "\n";
        for(auto& I : t.inds()) s << I << "\n";
        }
    else
        {
        s << t.inds();
        }
    if(not hasQNs(t)) s << "\n";
    if(not t.store()) 
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
matrixITensor(Matrix&& M, Index const& i1, Index const& i2)
    {
    auto res = ITensor({i1,i2},DenseReal{std::move(M.storage())});
    M.clear();
    return res;
    }

//Deprecated
ITensor
matrixTensor(Matrix&& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(Matrix,Index,Index) is deprecated in favor of matrixITensor(Matrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }

ITensor
matrixITensor(Matrix const& M, Index const& i1, Index const& i2)
    {
    return matrixITensor(Matrix(M),i1,i2);
    }

//Deprecated
ITensor
matrixTensor(Matrix const& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(Matrix,Index,Index) is deprecated in favor of matrixITensor(Matrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }

ITensor
matrixITensor(CMatrix&& M, Index const& i1, Index const& i2)
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

//Deprecated
ITensor
matrixTensor(CMatrix&& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(CMatrix,Index,Index) is deprecated in favor of matrixITensor(CMatrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }

ITensor
matrixITensor(CMatrix const& M, Index const& i1, Index const& i2)
    {
    return matrixITensor(CMatrix(M),i1,i2);
    }

//Deprecated
ITensor
matrixTensor(CMatrix const& M, Index const& i1, Index const& i2)
    {
    Global::warnDeprecated("matrixTensor(CMatrix,Index,Index) is deprecated in favor of matrixITensor(CMatrix,Index,Index)");
    return matrixITensor(M,i1,i2);
    }


ITensor
combiner(IndexSet const& inds, Args const& args)
    {
    auto itagset = getTagSet(args,"Tags","Link,CMB");

    if(not hasQNs(inds))
        {
        if(inds.empty()) Error("No indices passed to combiner");
        long rm = 1;
        for(const auto& i : inds) rm *= dim(i);
        //create combined index
        auto cind = Index(rm,itagset);
        //create new IndexSet with combined index in front
        auto newind = IndexSetBuilder(1+inds.r());
        newind.nextIndex(std::move(cind));
        for(auto& I : inds)
            newind.nextIndex(std::move(I));
        return ITensor(newind.build(),Combiner());
        }
    else if(hasQNs(inds))
        {
        if(inds.empty()) Error("No indices passed to combiner");

        auto cdir = Out;
        if(args.defined("IndexDir"))
            {
            cdir = toArrow(args.getInt("IndexDir"));
            }
        else
            {
            //If not specified by user, make combined Index
            //point Out unless all combined indices have In arrows.
            auto allin = true;
            for(auto& i : inds) 
                if(i.dir() != In)
                    {
                    allin = false;
                    break;
                    }
            if(allin) cdir = In;
            }


        auto C = QCombiner(inds);

        //Build the combined Index,
        //saving information about
        //how we're doing it in C
        struct QNm { QN q; long m = 1; };
        auto qms = vector<QNm>{};

        //Loop over all possible QN sectors that
        //can be formed by the combined indices
        for(auto I : C.range())
            {
            QNm qm;
            //For this sector, figure out the total QN (qm.q)
            //and combined sector size (qm.m)
            for(auto j : range(inds))
                {
                qm.q += inds[j].qn(1+I[j]) * inds[j].dir() * cdir;
                qm.m *= inds[j].blocksize0(I[j]);
                }

            size_t n = 0;
            for(; n < qms.size(); ++n) if(qms[n].q==qm.q) break;

            if(n < qms.size())
                {
                C.setBlockRange(I,n,qms[n].m,qm.m);
                qms[n].m += qm.m;
                }
            else 
                {
                C.setBlockRange(I,n,0,qm.m);
                qms.push_back(qm);
                }
            }

        auto cstore = stdx::reserve_vector<QNInt>(qms.size());
        for(auto n : range(qms)) 
            {
            cstore.emplace_back(qms[n].q,qms[n].m);
            }
        auto cind = Index(std::move(cstore),cdir,itagset);

        auto newind = IndexSetBuilder(1+inds.size());
        newind.nextIndex(std::move(cind));
        for(auto& I : inds) 
            {
            newind.nextIndex(dag(I));
            }
        return ITensor(newind.build(),std::move(C));
        }
    return ITensor();
    }

ITensor
combiner(std::vector<Index> const& inds,
         Args const& args)
    {
    return combiner(IndexSet(inds),args);
    }

ITensor
combiner(std::initializer_list<Index> inds,
         Args const& args)
    {
    return combiner(IndexSet(inds),args);
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
randomITensor(IndexSet const& inds)
    {
    return random(ITensor{inds});
    }

ITensor
randomITensor(QN q, IndexSet const& is, Args const& args)
    {
#ifdef DEBUG
    if(not hasQNs(is)) 
        Error("Cannot use randomITensor(QN,...) to create non-QN-conserving ITensor");
#endif
    ITensor T;
    if(args.getBool("Complex",false))
        {
        auto dat = QDenseCplx{is,q};
        T = ITensor(std::move(is),std::move(dat));
        T.generate(detail::quickranCplx);
        }
    else
        {
        auto dat = QDenseReal{is,q};
        T = ITensor(std::move(is),std::move(dat));
        T.generate(detail::quickran);
        }
    return T;
    }

QN
div(ITensor const& T) 
    { 
    if(not hasQNs(T)) Error("div(ITensor) not defined for non QN conserving ITensor");
    if(!T) Error("div(ITensor) not defined for unallocated IQTensor");
    return doTask(CalcDiv{T.inds()},T.store());
    }

QN
flux(ITensor const& T) 
    {
    return div(T);
    }

namespace detail {

IndexSet
moveToFront(IndexSet const& isf, IndexSet const& is)
    {
    auto rf = isf.r();
    auto r = is.r();

    if(rf >= r)
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n",isf," '...'\n");
        println("---------------------------------------------");
        Error(format("Wrong number of indices passed to permute (expected < %d, got %d)",r,rf));
        }

    auto iso = IndexSet(r);

    auto i = 0;
    for(auto& I : isf) 
        {
        if(!hasIndex(is,I))
            {
            println("---------------------------------------------");
            println("Tensor indices = \n",is,"\n");
            println("---------------------------------------------");
            println("Indices provided = \n",isf," '...'\n");
            println("---------------------------------------------");
            Error(format("Bad index passed to permute"));
            }
        iso[i] = I;
        i++;
        }

    auto j = rf;
    for(auto& J : is)
        {
        if(!hasIndex(isf,J))
            {
            iso[j] = J;
            j++;
            }
        }

    return iso;
    }

IndexSet 
moveToBack(IndexSet const& isb, IndexSet const& is)
    {
    auto rb = isb.r();
    auto r = is.r();

    if(rb >= r)
        {
        println("---------------------------------------------");
        println("Tensor indices = \n",is,"\n");
        println("---------------------------------------------");
        println("Indices provided = \n'...' ",isb,"\n");
        println("---------------------------------------------");
        Error(format("Wrong number of indices passed to permute (expected < %d, got %d)",r,rb));
        }

    auto iso = IndexSet(r);

    auto i = r-rb;
    for(auto& I : isb) 
        {
        if(!hasIndex(is,I))
            {
            println("---------------------------------------------");
            println("Tensor indices = \n",is,"\n");
            println("---------------------------------------------");
            println("Indices provided = \n'...' ",isb,"\n");
            println("---------------------------------------------");
            Error(format("Bad index passed to permute"));
            }
        iso[i] = I;
        i++;
        }

    auto j = 0;
    for(auto& J : is)
        {
        if(!hasIndex(isb,J))
            {
            iso[j] = J;
            j++;
            }
        }

    return iso;
    }

} //namespace detail

} //namespace itensor
