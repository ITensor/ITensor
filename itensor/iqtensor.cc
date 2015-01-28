//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"
#include "detail/printing.h"
#include "lapack_wrap.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::find;
using std::pair;
using std::make_pair;
using std::make_shared;
using std::move;

void
inverseBlockInd(long I,
                const IQIndexSet& is,
                vector<long>& ind)
    {
    auto r = int(ind.size());
    assert(r == is.r());
    for(int j = 0; j < r-1; ++j)
        {
        ind[j] = I % is[j].nindex();
        I = (I-ind[j])/is[j].nindex();
        }
    ind[r-1] = I;
    }

template<typename Indexable>
class IndexDim
    {
    const IQIndexSet& is_;
    const Indexable& ind_;
    public:

    IndexDim(const IQIndexSet& is,
             const Indexable& ind)
        :
        is_(is),
        ind_(ind)
        { }

    long
    size() const { return is_.r(); }

    long
    operator[](long j) const { return (is_[j])[ind_[j]].m(); }
    };

template<typename Indexable>
IndexDim<Indexable>
make_indexdim(const IQIndexSet& is, const Indexable& ind) 
    { return IndexDim<Indexable>(is,ind); }


IQTensor::
IQTensor(Complex val) 
    :
    scale_(1.)
    { 
    if(val.imag()==0)
        store_ = make_shared<ITDiag<Real>>(val.real());
    else
        store_ = make_shared<ITDiag<Complex>>(val);
    }

IQTensor::
IQTensor(const QN& q, vector<IQIndex>&& iqinds) 
	: 
    is_(move(iqinds)),
    store_(make_shared<IQTData<Real>>(is_,q)),
    div_(q),
    scale_(1.)
	{ }

class IQPlusEQ
    {
    Real fac_;
    const Permutation *P_ = nullptr;
    const IQIndexSet *is1_ = nullptr,
                     *is2_ = nullptr;
    public:

    IQPlusEQ(Real fac) : fac_(fac) { }

    IQPlusEQ(const Permutation& P,
             const IQIndexSet& is1,
             const IQIndexSet& is2,
             Real fac)
        :
        fac_(fac),
        P_(&P),
        is1_(&is1),
        is2_(&is2)
        { }

    ITResult
    operator()(IQTData<Real>& a1,
               const IQTData<Real>& a2);

    };

ITResult IQPlusEQ::
operator()(IQTData<Real>& A,
           const IQTData<Real>& B)
    {
#ifdef DEBUG
    if(A.data.size() != B.data.size()) Error("Mismatched sizes in plusEq");
#endif
    if(!P_)
        {
        LAPACK_INT inc = 1;
        LAPACK_INT size = A.data.size();
        daxpy_wrapper(&size,&fac_,B.data.data(),&inc,A.data.data(),&inc);
        }
    else
        {
        auto r = is1_->r();
        vector<long> Ablock(r,0),
                     Bblock(r,0);
        Range Arange,
              Brange;
        for(const auto& aio : A.offsets)
            {
            inverseBlockInd(aio.block,*is1_,Ablock);
            for(int i = 0; i < r; ++i)
                Bblock[i] = Ablock[P_->dest(i)];
            Arange.init(make_indexdim(*is1_,Ablock));
            Brange.init(make_indexdim(*is2_,Bblock));
            const auto* bblock = B.getBlock(*is2_,Bblock);

            auto aref = make_tensorref(A.data.data()+aio.offset,Arange),
                 bref = make_tensorref(bblock,Brange);
            auto f = fac_;
            auto add = [f](Real& r1, Real r2) { r1 += f*r2; };
            reshape(bref,*P_,aref,add);
            }
        }
    return ITResult();
    }

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    if(!*this) Error("Calling += on default constructed IQTensor");
    if(!other) Error("Right-hand-side of IQTensor += is default constructed");
    if(this == &other) return operator*=(2.);
    if(this->scale_.isZero()) return operator=(other);

    Permutation P(is_.size());
    try {
        detail::calc_permutation(other.is_,is_,P);
        }
    catch(const ITError& e)
        {
        Print(*this);
        Print(other);
        Error("IQTensor::operator+=: different IQIndex structure");
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

    if(isTrivial(P))
        {
        applyFunc<IQPlusEQ>(store_,other.store_,{scalefac});
        }
    else
        {
        applyFunc<IQPlusEQ>(store_,other.store_,{P,is_,other.is_,scalefac});
        }


    return *this;
    }

IQTensor& IQTensor::
operator*=(Real fac)
    {
    scale_ *= fac;
    return *this;
    }

IQTensor& IQTensor::
operator/=(Real fac)
    {
    scale_ /= fac;
    return *this;
    }

IQTensor& IQTensor::
operator*=(const LogNumber& lgnum)
    {
    scale_ *= lgnum;
    return *this;
    }

IQTensor& IQTensor::
operator-=(const IQTensor& o)
    { 
    if(this == &o) { operator*=(0); return *this; }
    IQTensor oth(o);
    oth *= -1;
    return operator+=(oth);
    }

class QContract
    {
    const Label &Aind_,
                &Bind_;

    const IQIndexSet &Ais_,
                     &Bis_;
    QN Cdiv_;
    IQIndexSet Cis_;
    Real scalefac_ = -1;
    public:

    QContract(const IQIndexSet& Ais,
              const Label& Aind,
              const IQIndexSet& Bis,
              const Label& Bind,
              const QN& Cdiv)
        :
        Aind_(Aind),
        Bind_(Bind),
        Ais_(Ais),
        Bis_(Bis),
        Cdiv_(Cdiv)
        { }

    IQIndexSet
    newIndexSet() { return move(Cis_); }

    Real
    scalefac() const { return scalefac_; }

    template<typename T>
    ITResult
    operator()(const IQTData<T>& d1,
               const IQTData<T>& d2);

    }; //QContract



template<typename T>
ITResult QContract::
operator()(const IQTData<T>& A,
           const IQTData<T>& B)
    {
    //compute C index set and sort:
    contractIS(Ais_,Aind_,Bis_,Bind_,Cis_,true);
    auto res = make_newdata<IQTData<Real>>(Cis_,Cdiv_);
    auto& C = *res;
    //println("Cis_ = \n",Cis_);

    auto rA = Ais_.r(),
         rB = Bis_.r(),
         rC = Cis_.r();

    Label AtoB(rA,-1),
          AtoC(rA,-1),
          BtoC(rB,-1);
    Label Cind(rC,0);
    for(size_t ic = 0; ic < rC; ++ic)
        {
        auto j = findindex(Ais_,Cis_[ic]);
        if(j >= 0)
            {
            Cind[ic] = Aind_[j];
            AtoC[j] = ic;
            }
        else
            {
            j = findindex(Bis_,Cis_[ic]);
            Cind[ic] = Bind_[j];
            BtoC[j] = ic;
            }
        }
    for(int ia = 0; ia < rA; ++ia)
    for(int ib = 0; ib < rB; ++ib)
        if(Aind_[ia] == Bind_[ib])
            {
            AtoB[ia] = ib;
            break;
            }
    
    detail::GCounter couB(rB);
    vector<long> Ablock(rA,0),
                 Cblock(rC,0);
    Range Arange,
          Brange,
          Crange;
    for(const auto& aio : A.offsets)
        {
        inverseBlockInd(aio.block,Ais_,Ablock);
        //PRI(Ablock);
        couB.reset();
        for(int ib = 0; ib < rB; ++ib)
            couB.setInd(ib,0,Bis_[ib].nindex()-1);
        for(int iA = 0; iA < rA; ++iA)
            {
            auto ival = Ablock[iA];
            if(AtoB[iA] != -1) couB.setInd(AtoB[iA],ival,ival);
            if(AtoC[iA] != -1) Cblock[AtoC[iA]] = ival;
            }
        for(;couB.notDone(); ++couB)
            {
            auto* bblock = B.getBlock(Bis_,couB.i);
            if(!bblock) continue;

            //Finish making Cblock
            for(int ib = 0; ib < rB; ++ib)
                if(BtoC[ib] != -1) Cblock[BtoC[ib]] = couB.i[ib];

            auto* cblock = C.getBlock(Cis_,Cblock);
            assert(cblock != nullptr);

            //PRI(couB.i);
            //println("aoff = ",aio.offset);
            //println("boff = ",boff);
            //println("coff = ",coff);

            Arange.init(make_indexdim(Ais_,Ablock));
            Brange.init(make_indexdim(Bis_,couB.i));
            Crange.init(make_indexdim(Cis_,Cblock));

            auto aref = make_tensorref(A.data.data()+aio.offset,Arange),
                 bref = make_tensorref(bblock,Brange);
            auto cref= make_tensorref(cblock,Crange);

            //PRI(Aind_);
            //PRI(Bind_);
            //PRI(Cind);
            //Print(Arange);
            //Print(Brange);
            //Print(Crange);
            //println("----------Calling contractloop--------------");
            //Print(aref);
            //Print(bref);
            contract(aref,Aind_,bref,Bind_,cref,Cind);
            //Print(cref);

            } //for couB
        } //for A.offsets

    scalefac_ = 0;
    for(auto elt : C.data) scalefac_ += elt*elt;
    scalefac_ = std::sqrt(scalefac_);
    if(scalefac_ != 0)
        {
        for(auto& elt : C.data) elt /= scalefac_;
        }

    return move(res);
    }

IQTensor& IQTensor::
operator*=(const IQTensor& other)
    {
    if(!(*this) || !other)
        Error("Default constructed IQTensor in product");

    if(this == &other)
        return operator=( IQTensor(sqr(norm(*this))) );

    const auto& Lis = is_;
    const auto& Ris = other.is_;

    auto checkDirs = 
    [&Lis,&Ris](const IQIndex& li, const IQIndex& ri)
        {
        if(li.dir() == ri.dir())
            {
            println("-------------------------");
            println("Left indices = \n",Lis);
            println("-------------------------");
            println("Right indices = \n",Ris);
            println("-------------------------");
            println("IQIndex from left = \n",li);
            println("IQIndex from right = \n",ri);
            Error("Incompatible arrow directions in IQTensor contraction");
            }
        };
    Label Lind,
          Rind;
    computeLabels(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind,checkDirs);

    //PRI(Lind);
    //PRI(Rind);

    auto C = applyFunc<QContract>(store_,other.store_,{Lis,Lind,Ris,Rind,div_+other.div_});

    is_ = C.newIndexSet();

    div_ += other.div_;

    scale_ *= other.scale_;
    if(C.scalefac() > 0) scale_ *= C.scalefac();

    return *this;
    }

IQTensor& IQTensor::
dag()
    {
    if(isComplex(*this))
        {
        Error("Not implemented");
        }
    is_.dag();
    div_ = -div_;
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
    operator()(IQTData<T>& d) const
        {
        //use BLAS algorithm?
        for(auto& elt : d.data)
            elt *= r_;
        return ITResult();
        }
    };

void IQTensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    applyFunc<MultReal>(store_,{scale_.real0()});
    scale_ = newscale;
    }

class ToITensor
    {
    ITensor res;
    const IQIndexSet& is_;
    const LogNumber& scale_;
    public:

    ToITensor(const IQIndexSet& is,
              const LogNumber& scale)
        :
        is_(is),
        scale_(scale)
        { }

    explicit
    operator ITensor() { return std::move(res); }

    template<typename T>
    ITResult
    operator()(const IQTData<T>& d)
        {
        auto r = is_.r();
        auto nd = make_newdata<ITDense<T>>(area(is_),0);
        auto *pd = d.data.data();
        auto *pn = nd->data.data();
        vector<long> block(r,0);
        detail::GCounter C(r);
        for(const auto& io : d.offsets)
            {
            inverseBlockInd(io.block,is_,block);
            for(long j = 0; j < r; ++j)
                {
                long start = 0;
                for(long b = 0; b < block[j]; ++b)
                    start += is_[j][b].m();
                C.setInd(j,start,start+is_[j][block[j]].m()-1);
                }
            for(; C.notDone(); ++C)
                {
                pn[ind(is_,C.i)] = pd[io.offset+C.ind];
                }
            }
        vector<Index> inds(r);
        for(long j = 0; j < r; ++j) inds[j] = is_[j];
        res = ITensor(IndexSet(std::move(inds)),std::move(nd),scale_);
        return ITResult();
        }
    };

ITensor
toITensor(const IQTensor& T)
    {
    return ITensor(applyFunc<ToITensor>(T.data(),{T.inds(),T.scale()}));
    }

struct IsComplex
    {
    bool res = false;
    IsComplex() { }

    operator bool() const { return res; }

    ITResult
    operator()(const IQTData<Real>& d)
        {
        res = false;
        return ITResult();
        }

    ITResult
    operator()(const IQTData<Complex>& d)
        {
        res = true;
        return ITResult();
        }
    };


bool
isComplex(const IQTensor& T)
    {
    return applyFunc<IsComplex>(T.data());
    }

IQIndex
findIQInd(const IQTensor& T, const Index& i)
    {
    for(const IQIndex& J : T.indices())
        if(hasindex(J,i)) return J;
    Print(T.indices());
    Print(i);
    throw ITError("Index i not found in any of T's IQIndices");
    return IQIndex();
    }


Arrow
dir(const IQTensor& T, const IQIndex& I)
	{
    for(const IQIndex& J : T.indices())
        if(I == J) return J.dir();
    Error("dir: IQIndex not found");
    return Out;
	}

class NormNoScale
    {
    Real nrm_;
    public:

    NormNoScale() : nrm_(0) { }

    operator Real() const { return nrm_; }

    template<typename T>
    ITResult
    operator()(const IQTData<T>& d) { calc(d); return ITResult(); }

    private:

    template<typename T>
    void
    calc(const T& d)
        {
        for(const auto& elt : d.data)
            {
            nrm_ += std::norm(elt);
            }
        nrm_ = std::sqrt(nrm_);
        }
    };

Real
norm(const IQTensor& T)
    {
#ifdef DEBUG
    if(!T) Error("ITensor is default initialized");
#endif
    return T.scale().real0() *
           applyFunc<NormNoScale>(T.data());
    }

IQTensor
randomize(IQTensor T, const Args& args)
    {
    T.generate(detail::quickran);
    return T;
    }

bool
isZero(const IQTensor& T, const Args& args)
    {
    Error("Not implemented");
    //if(T.empty()) return true;
    ////done with all fast checks
    //if(args.getBool("Fast",false)) return false;
    //for(const ITensor& t : T.blocks())
    //    {
    //    if(!isZero(t)) return false;
    //    }
    return true;
    }

struct PrintIQT
    {
    std::ostream& s_;
    const LogNumber& x_;
    const IQIndexSet& is_;

    PrintIQT(std::ostream& s,
             const LogNumber& x,
             const IQIndexSet& is)
        : s_(s), x_(x), is_(is)
        { }

    template<typename T>
    ITResult
    operator()(const IQTData<T>& d) const;
    };

template<typename T>
ITResult PrintIQT::
operator()(const IQTData<T>& d) const
    {
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "(omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) return ITResult();
        
    vector<long> block(rank,0);
    auto blockIndex = [&block,this](long i)->const Index& { return (this->is_[i])[block[i]]; };

    Range brange;
    detail::GCounter C(rank);
    for(const auto& io : d.offsets)
        {
        //Determine block indices (where in the IQIndex space
        //this non-zero block is located)
        inverseBlockInd(io.block,is_,block);
        //Print Indices of this block
        for(size_t i = 0; i < rank; ++i)
            {
            s_ << blockIndex(i) << "<" << is_[i].dir() << "> ";
            }
        s_ << "\n";
        //Wire up GCounter with appropriate dims
        C.reset();
        for(int i = 0; i < rank; ++i)
            C.setInd(i,0,blockIndex(i).m()-1);
        for(auto os = io.offset; C.notDone(); ++C, ++os)
            {
            auto val = scalefac*d.data[os];
            if(std::norm(val) > Global::printScale())
                {
                s_ << "  (";
                for(auto ii = C.i.mini(); ii <= C.i.maxi(); ++ii)
                    {
                    s_ << (1+C.i(ii));
                    if(ii < C.i.maxi()) s_ << ",";
                    }
                s_ << ") ";
                detail::printVal(s_,val);
                }
            }
        }

    return ITResult();
    }

std::ostream&
operator<<(std::ostream& s, const IQTensor& T)
    {
	s << "/--------------IQTensor--------------\n";
    if(T)
        {
        s << "r=" << T.r() << " div=" << div(T) << " log(scale)=" << T.scale().logNum() << "\n";
        s << T.inds();
        //Checking whether std::ios::floatfield is set enables 
        //printing the contents of an ITensor when using the printf
        //format string %f (or another float-related format string)
        const bool ff_set = (std::ios::floatfield & s.flags()) != 0;

        if(ff_set || Global::printdat())
            {
            s << "\n";
            applyFunc<PrintIQT>(T.data(),{s,T.scale(),T.inds()});
            }
        }
    else
        {
        s << "(default constructed)\n";
        }
	s << "\\------------------------------------\n\n";
    return s;
    }


}; //namespace itensor
