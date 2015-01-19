//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"

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
inverseInd(long I,
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
        inverseInd(aio.ind,Ais_,Ablock);
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
            auto boff = B.getOffset(couB.i,[this](long i){ return this->Bis_[i].nindex(); });
            if(boff < 0) continue;

            //Finish making Cblock
            for(int ib = 0; ib < rB; ++ib)
                if(BtoC[ib] != -1) Cblock[BtoC[ib]] = couB.i.fast(ib);

            auto coff = C.getOffset(Cblock,[this](long i){ return this->Cis_[i].nindex(); });
            assert(coff != -1);

            PRI(couB.i);
            println("aoff = ",aio.offset);
            println("boff = ",boff);
            println("coff = ",coff);
            println();

            Arange.init(make_indexdim(Ais_,Ablock));
            Brange.init(make_indexdim(Bis_,couB.i));
            Crange.init(make_indexdim(Cis_,Cblock));

            auto aref = make_tensorref(A.data.data()+aio.offset,Arange),
                 bref = make_tensorref(B.data.data()+boff,Brange);
            auto cref= make_tensorref(C.data.data()+coff,Crange);

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
    Error("Not implemented");
    div_ = -div_;
    return *this;
    }

ITensor
toITensor(const IQTensor& T)
    {
    Error("toITensor not implemented");
    return ITensor();
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


}; //namespace itensor
