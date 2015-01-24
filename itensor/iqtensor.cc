//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"
#include "detail/printing.h"

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

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    Error("Not implemented");

    //
    // Idea to implement:
    // - Compute permutation P from is_ to other.is_
    // - If P trivial, do a straight daxpy on the data
    // - If P not trivial, loop over blocks of *this
    //   and get pointer to corresponding block of other 
    //   by applying P and add those blocks.
    //   (Instead of looping over a GCounter
    //   to visit all non-zero blocks, may be faster
    //   to iterate through offset_ and calculate block indices
    //   from position of non-negative offset_ elements.)
    //

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
    
    detail::GCounter couB(0,rB-1,0);
    vector<long> Ablock(rA,0),
                 Cblock(rC,0);
    Range Arange,
          Brange,
          Crange;
    for(const auto& aio : A.offsets)
        {
        inverseBlockInd(aio.block,Ais_,Ablock);
        PRI(Ablock);
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

            PRI(couB.i);
            //println("aoff = ",aio.offset);
            //println("boff = ",boff);
            //println("coff = ",coff);
            println();

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
            //Print(cref);
            contractloop(aref,Aind_,bref,Bind_,cref,Cind);
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

    PRI(Lind);
    PRI(Rind);

    auto C = applyFunc<QContract>(store_,other.store_,{Lis,Lind,Ris,Rind,div_+other.div_});

    is_ = C.newIndexSet();

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
        solo();
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

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("IQTensor MultReal not implemented for ITData type."); return ITResult(); }
    };

void IQTensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    solo();
    scale_ /= newscale;
    applyFunc<MultReal>(store_,{scale_.real0()});
    scale_ = newscale;
    }

void IQTensor::
solo()
	{
    if(!store_.unique()) store_ = store_->clone();
    }

ITensor
toITensor(const IQTensor& T)
    {
    Error("toITensor not implemented");
    return ITensor();
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

    template<typename T>
    ITResult
    operator()(const T& d) const { Error("Function not implemented."); return ITResult(); }
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
    detail::GCounter C(0,rank-1,0);
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
    s << "r=" << T.r() << ", log(scale)[incl in elems]=" << T.scale().logNum() << "\n";
    s << T.inds();

    //Checking whether std::ios::floatfield is set enables 
    //printing the contents of an ITensor when using the printf
    //format string %f (or another float-related format string)
    const bool ff_set = (std::ios::floatfield & s.flags()) != 0;

    if(ff_set || Global::printdat())
        {
        if(T) applyFunc<PrintIQT>(T.data(),{s,T.scale(),T.inds()});
        else           s << " (default constructed)}\n";
        }
	s << "\\------------------------------------\n\n";
    return s;
    }


}; //namespace itensor
