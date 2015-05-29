//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/matrix/lapack_wrap.h"
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
  : is_(i1),
    store_(std::make_shared<ITReal>(i1.m(),0.)),
    scale_(1.)
	{ }


template<>
ITensor::
ITensorT(const Index& i1,const Index& i2) 
  : is_(i1,i2),
    store_(std::make_shared<ITReal>(i1.m()*i2.m(),0.)),
    scale_(1.)
	{ }
    
template<>
ITensor::
ITensorT(Cplx val) 
  : scale_(1.)
    { 
    if(val.imag() == 0)
        store_ = std::make_shared<ITReal>(1,val.real());
    else
        store_ = std::make_shared<ITCplx>(1,val);
    //if(val.imag() == 0)
    //    store_ = std::make_shared<ITDiag<Real>>(val.real());
    //else
    //    store_ = std::make_shared<ITDiag<Cplx>>(val);
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
    computeLabels(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind);

    auto nstore = A.store();
    auto C = doTask(Contract<Index>{Lis,Lind,Ris,Rind},nstore,B.store());

    auto nscale = A.scale() * B.scale();
    if(!std::isnan(C.scalefac)) nscale *= C.scalefac;

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    A = ITensor(C.Nis,std::move(nstore),nscale);

    return A;
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

ITensor&
operator+=(ITensor& A, const ITensor& B)
    {
    if(!A) Error("Calling += on default constructed ITensor");
    if(!B) Error("Right-hand-side of ITensor += is default constructed");
    if(&A == &B) return operator*=(A,2.);
    if(A.scale().isZero()) return A.operator=(B);

    PlusEQ<Index>::permutation P(A.inds().size());
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
        doTask(PlusEQ<Index>{scalefac},A.store(),B.store());
        }
    else
        {
        doTask(PlusEQ<Index>{P,A.inds(),B.inds(),scalefac},A.store(),B.store());
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

//template<>
//void ITensor::
//scaleOutNorm()
//    {
//    auto nrm = doTask<Real>(NormNoScale<Index>{is_},store_);
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

template<>
ITensor& ITensor::
conj()
    {
    doTask(Conj{},store_);
    return *this;
    }

template<>
ITensor& ITensor::
takeReal()
    {
    doTask(TakeReal{},store_);
    return *this;
    }

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
        doTask(PrintIT<Index>{s,t.scale(),t.inds(),print_data},t.store());
        }
    return s;
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
           doTask<Real>(NormNoScale<Index>{T.inds()},T.store());
    }

ITensor
conj(ITensor T)
    {
    T.conj();
    return T;
    }

bool
isComplex(const ITensor& t)
    {
    return doTask<bool>(CheckComplex{},t.store());
    }

Cplx
sumelsC(const ITensor& t)
    {
    auto z = doTask<Cplx>(SumEls<Index>{t.inds()},t.store());
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
    auto type = StorageType::Null;
    s.read((char*)&type,sizeof(type));
    ITensor::storage_ptr p;
    if(type==StorageType::Null) { /*intentionally left blank*/  }
    else if(type==StorageType::ITReal) { p = readType<ITReal>(s); }
    else if(type==StorageType::ITCplx) { p = readType<ITCplx>(s); }
    else if(type==StorageType::ITCombiner) { p = readType<ITCombiner>(s); }
    else if(type==StorageType::ITDiagReal) { p = readType<ITDiag<Real>>(s); }
    else if(type==StorageType::ITDiagCplx) { p = readType<ITDiag<Cplx>>(s); }
    else
        {
        Error("Unrecognized type when reading ITensor from istream");
        }
    T = ITensor(std::move(is),std::move(p),scale);
    }

void
write(std::ostream& s, const ITensor& T)
    {
    write(s,T.inds());
    write(s,T.scale());
    if(T) 
        doTask(Write{s},T.store());
    else 
        write(s,StorageType::Null);
    }

} //namespace itensor
