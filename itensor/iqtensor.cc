//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"
#include "detail/printing.h"
#include "matrix/lapack_wrap.h"
#include "tensor/contract.h"
#include "count.h"

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
using std::shared_ptr;
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

IQTensor::
IQTensor(const QN& q,
         IQIndexSet&& iset,
         NewData nd,
         LogNumber scale)
    :
    is_(move(iset)),
    store_(move(nd)),
    div_(q),
    scale_(scale)
    {
    }

class IQPlusEQ : public RegisterFunc<IQPlusEQ>
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

    void
    operator()(IQTData<Real>& a1,
               const IQTData<Real>& a2);

    };

void IQPlusEQ::
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

            auto aref = makeTensorRef(A.data.data()+aio.offset,Arange);
            auto bref = makeTensorRef(bblock,Brange);
            auto f = fac_;
            auto add = [f](Real& r1, Real r2) { r1 += f*r2; };
            permute(bref,*P_,aref,add);
            }
        }
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
        calc_permutation(other.is_,is_,P);
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
        applyFunc<IQPlusEQ>(store_,other.store_,scalefac);
        }
    else
        {
        applyFunc<IQPlusEQ>(store_,other.store_,P,is_,other.is_,scalefac);
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

template<typename T>
void
permuteIQ(const Permutation& P,
          const IQIndexSet& Ais,
          const QN& div,
          const IQTData<T>& dA,
          IQIndexSet& Bis,
          shared_ptr<IQTData<T>>& pdB)

    {
#ifdef DEBUG
    if(isTrivial(P)) Error("Calling permuteIQ for trivial Permutation");
#endif
    auto r = Ais.r();
    vector<IQIndex> bind(r);
    for(auto i : count(r))
        {
        bind.at(P.dest(i)) = Ais[i];
        }
    Bis = IQIndexSet(std::move(bind));
    pdB = make_shared<IQTData<T>>(Bis,div);

    vector<long> Ablock(r,-1),
                 Bblock(r,-1);
    Range Arange,
          Brange;
    for(auto aio : dA.offsets)
        {
        //Compute bi, new block index of blk
        inverseBlockInd(aio.block,Ais,Ablock);
        for(auto j : index(Ablock)) 
            Bblock.at(P.dest(j)) = Ablock[j];
        Arange.init(make_indexdim(Ais,Ablock));
        Brange.init(make_indexdim(Bis,Bblock));

        auto* bblock = pdB->getBlock(Bis,Bblock);
        auto aref = makeTensorRef(dA.data.data()+aio.offset,Arange);
        auto bref = makeTensorRef(bblock,Brange);
        permute(aref,P,bref);
        }
    }

class QContract : public RegisterFunc<QContract>
    {
    const Label &Aind_,
                &Bind_;

    const IQIndexSet &Ais_,
                     &Bis_;
    QN Cdiv_;
    IQIndexSet Nis_;
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
    newIndexSet() { return move(Nis_); }

    Real
    scalefac() const { return scalefac_; }

    template<typename T>
    void
    operator()(const IQTData<T>& d1,
               const IQTData<T>& d2);

    void
    operator()(const IQTData<Real>& d,
               const ITCombiner& C)
        {
        combine(d,Ais_,Bis_,true);
        }
    void
    operator()(const ITCombiner& C,
               const IQTData<Real>& d)
        { 
        combine(d,Bis_,Ais_,false);
        }

    private:

    void
    combine(const IQTData<Real>& d,
            const IQIndexSet& dis,
            const IQIndexSet& Cis,
            bool own_data);

    }; //QContract

void QContract::
combine(const IQTData<Real>& d,
        const IQIndexSet& dis,
        const IQIndexSet& Cis,
        bool own_data)
    {
    //cind is "combined index"
    const auto& cind = Cis[0];
    //check if d has combined index i.e. we are "uncombining"
    auto jc = findindex(dis,cind);

    Permutation P;
    if(jc > 0) //we are uncombining, but cind not at front
        {
        P = Permutation(dis.r());
        P.setFromTo(jc,0);
        long ni = 1;
        for(auto j : count(dis.r()))
            if(j != jc) P.setFromTo(j,ni++);
        jc = 0;
        }
    else if(jc < 0)
        {
        //check locations of Cis[1], Cis[2], ...
        //Check if Cis[1],Cis[2],... are grouped together (contiguous)
        //all at front, and in same order as on combiner
        bool front_contig = true;
        for(auto j = 0, c = 1; c < Cis.r() && j < dis.r(); ++j,++c)
            if(dis[j] != Cis[c])
                {
                front_contig = false;
                break;
                }
        if(!front_contig) //if !front_contig, need to permute
            {
            P = Permutation(dis.r());
            //Set P destination values to -1 to mark
            //indices that need to be assigned destinations:
            for(auto i : index(P)) P.setFromTo(i,-1);

            //permute combined indices to the front, in same
            //order as in Cis:
            long ni = 0;
            for(auto c : count(1,Cis.r()))
                {
                auto j = findindex(dis,Cis[c]);
                if(j < 0) 
                    {
                    println("IQIndexSet of dense tensor =\n  ",dis);
                    println("IQIndexSet of combiner/delta =\n  ",Cis);
                    println("Missing IQIndex: ",Cis[c]);
                    println("jc = ",jc);
                    Error("IQCombiner: missing IQIndex");
                    }
                P.setFromTo(j,ni++);
                }
            for(auto j : index(P))
                {
                if(P.dest(j) == -1) P.setFromTo(j,ni++);
                }
            }
        }

    std::shared_ptr<IQTData<Real>> pnd;
    IQTData<Real>* nd = nullptr;
    if(P) 
        {
        permuteIQ(P,dis,Cdiv_,d,Nis_,pnd);
        nd = pnd.get();
        }
    auto& Pis = (P ? Nis_ : dis);

    if(jc == 0) //has cind at front, we are "uncombining"
        {
        auto newr = Pis.r()+Cis.r()-2;
        auto offset = Cis.r()-1;
        vector<IQIndex> newind(newr);
        for(auto j : count(offset))
            newind.at(j) = Cis[1+j];
        for(auto j : count(Pis.r()-1))
            newind.at(offset+j) = Pis[1+j];
        Nis_ = IQIndexSet(move(newind));
        }
    else //we are "combining"
        {
        auto newr = Pis.r()-Cis.r()+2;
        auto offset = Cis.r()-2;
        vector<IQIndex> newind(newr);
        newind.front() = cind;
        for(auto j : count(1,newr))
            newind.at(j) = Pis[offset+j];
        Nis_ = IQIndexSet(move(newind));
        }

    //Only need to modify d if Cis.r() > 2.
    //If Cis.r()==2 just swapping one index for another
    if(Cis.r() > 2)
        {
        if(!nd)
            {
            if(own_data) 
                {
                nd = &modifyData(d);
                }
            else
                {
                pnd = std::make_shared<IQTData<Real>>(d);
                nd = pnd.get();
                }
            }
        nd->updateOffsets(Nis_,Cdiv_);
        }

    if(pnd) setNewData(std::move(pnd));
    }


template<typename T>
void QContract::
operator()(const IQTData<T>& A,
           const IQTData<T>& B)
    {
    //compute new index set (Nis_):
    contractIS(Ais_,Aind_,Bis_,Bind_,Nis_,true);

    //Allocate storage for C
    auto nd = makeNewData<IQTData<Real>>(Nis_,Cdiv_);
    auto& C = *nd;

    auto rA = Ais_.r(),
         rB = Bis_.r(),
         rC = Nis_.r();

    Label AtoB(rA,-1),
          AtoC(rA,-1),
          BtoC(rB,-1);
    Label Cind(rC,0);
    for(auto ic : count(rC))
        {
        auto j = findindex(Ais_,Nis_[ic]);
        if(j >= 0)
            {
            Cind[ic] = Aind_[j];
            AtoC[j] = ic;
            }
        else
            {
            j = findindex(Bis_,Nis_[ic]);
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
    //Loop over blocks of A (labeled by elements of A.offsets)
    for(const auto& aio : A.offsets)
        {
        //Reconstruct indices labeling this block of A, put into Ablock
        inverseBlockInd(aio.block,Ais_,Ablock);
        //Reset couB to run over indices of B (at first)
        couB.reset();
        for(int ib = 0; ib < rB; ++ib)
            couB.setInd(ib,0,Bis_[ib].nindex()-1);
        for(int iA = 0; iA < rA; ++iA)
            {
            auto ival = Ablock[iA];
            //Restrict couB to be fixed for indices of B contracted with A
            if(AtoB[iA] != -1) couB.setInd(AtoB[iA],ival,ival);
            //Begin computing elements of Cblock(=destination of this block-block contraction)
            if(AtoC[iA] != -1) Cblock[AtoC[iA]] = ival;
            }
        //Loop over blocks of B which contract with current block of A
        for(;couB.notDone(); ++couB)
            {
            //Check whether B contains non-zero block for this setting of couB
            //TODO: check whether block is present by computing its QN flux,
            //      should be faster than calling getBlock
            auto* bblock = B.getBlock(Bis_,couB.i);
            if(!bblock) continue;

            //Finish making Cblock index array
            for(int ib = 0; ib < rB; ++ib)
                if(BtoC[ib] != -1) Cblock[BtoC[ib]] = couB.i[ib];

            auto* cblock = C.getBlock(Nis_,Cblock);
            assert(cblock != nullptr);

            //Construct range objects for aref,bref,cref
            //using IndexDim helper objects
            Arange.init(make_indexdim(Ais_,Ablock));
            Brange.init(make_indexdim(Bis_,couB.i));
            Crange.init(make_indexdim(Nis_,Cblock));

            //"Wire up" TensorRef's pointing to blocks of A,B, and C
            //we are working with
            auto aref = makeTensorRef(A.data.data()+aio.offset,Arange),
                 bref = makeTensorRef(bblock,Brange);
            auto cref= makeTensorRef(cblock,Crange);

            //Compute aref*bref=cref
            contract(aref,Aind_,bref,Bind_,cref,Cind);

            } //for couB
        } //for A.offsets

    //Compute new scalefac_ from C.data
    scalefac_ = 0;
    for(auto elt : C.data) scalefac_ += elt*elt;
    scalefac_ = std::sqrt(scalefac_);
    //Rescale C by scalefac_
    if(scalefac_ != 0)
        {
        for(auto& elt : C.data) elt /= scalefac_;
        }
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

    auto qcres = 
    applyFunc<QContract>(store_,other.store_,Lis,Lind,Ris,Rind,div_+other.div_);

    is_ = qcres.newIndexSet();

    div_ += other.div_;

    scale_ *= other.scale_;
    if(qcres.scalefac() > 0) scale_ *= qcres.scalefac();

    return *this;
    }

struct AddITensor : RegisterFunc<AddITensor>
    {
    const IQIndexSet& iqis;
    const IndexSet& is;
    const vector<long>& block_ind;
    const Permutation& P;
    Real fac = 0;
    AddITensor(const IQIndexSet& iqis_,
               const IndexSet& is_,
               const vector<long>& block_ind_,
               const Permutation& P_,
               Real scalefac_)
        :
        iqis(iqis_),
        is(is_),
        block_ind(block_ind_),
        P(P_),
        fac(scalefac_)
        { }

    void
    operator()(IQTData<Real>& d, const ITReal& t)
        {
        Range drange;
        drange.init(make_indexdim(iqis,block_ind));
        auto* dblock = d.getBlock(iqis,block_ind);

        auto dref = makeTensorRef(dblock,drange);
        auto tref = makeTensorRef(t.data(),is);
        auto f = fac;
        auto add = [f](Real& r1, Real r2) { r1 += f*r2; };
        permute(tref,P,dref,add);
        }
    };

void
calcDiv(QN& d, const IQIndexSet& is, const vector<long>& block_ind)
    {
    d = QN();
    for(auto i : count(is.r())) { d += is[i].dir()*is[i].qn(1+block_ind[i]); }
    }

IQTensor& IQTensor::
operator+=(const ITensor& t)
    {
    if(!t) Error("IQTensor+=ITensor: r.h.s. ITensor is default constructed");
    if(!is_) Error("Calling IQTensor+= on default constructed ITensor");
    auto rank = r();
#ifdef DEBUG
    if(t.r() != rank) Error("Mismatched number of indices in IQTensor+=ITensor");
#endif

    Permutation P(rank);
    vector<long> block_ind(rank);
    for(auto i : count(rank))
    for(auto I : count(rank))
        {
        auto j = findindex(is_[I],t.inds()[i]);
        if(j > 0)
            {
            block_ind[I] = (j-1);
            P.setFromTo(i,I);
            break;
            }
        }

    if(!store_) 
        {
        //allocate data to add this ITensor into
        calcDiv(div_,is_,block_ind);
        if(!isComplex(t)) store_ = make_shared<IQTData<Real>>(is_,div_);
        else              Error("Initializing complex IQTensor in +=ITensor not yet implemented");
        }
#ifdef DEBUG
    else
        {
        QN q;
        calcDiv(q,is_,block_ind);
        if(q != div_) Error("IQTensor+=ITensor, ITensor has incompatible QN flux/divergence");
        }
#endif

    Real scalefac = 1;
    if(scale_.magnitudeLessThan(t.scale())) scaleTo(t.scale()); 
    else                                    scalefac = (t.scale()/scale_).real();

    applyFunc<AddITensor>(store_,t.data(),is_,t.inds(),block_ind,P,scalefac);
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

class MultReal : public RegisterFunc<MultReal>
    {
    Real r_;
    public:
    MultReal(Real r)
        : r_(r)
        { }

    template<typename T>
    void
    operator()(IQTData<T>& d) const
        {
        //use BLAS algorithm?
        for(auto& elt : d.data)
            elt *= r_;
        }
    };

void IQTensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an ITensor to a 0 scale");
    scale_ /= newscale;
    applyFunc<MultReal>(store_,scale_.real0());
    scale_ = newscale;
    }

class ToITensor : public RegisterFunc<ToITensor,ITensor>
    {
    const IQIndexSet& is_;
    const LogNumber& scale_;
    public:

    ToITensor(const IQIndexSet& is,
              const LogNumber& scale)
        :
        is_(is),
        scale_(scale)
        { }

    ITensor
    operator()(const IQTData<Real>& d)
        {
        auto r = is_.r();
        auto nd = ITReal(area(is_),0);
        auto *pd = d.data.data();
        auto *pn = nd.data();
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
        return ITensor(IndexSet(std::move(inds)),std::move(nd),scale_);
        }
    };

ITensor
toITensor(const IQTensor& T)
    {
    return applyFunc<ToITensor>(T.data(),T.inds(),T.scale());
    }

struct IsComplex : RegisterFunc<IsComplex,bool>
    {
    bool
    operator()(const IQTData<Complex>& d) { return true; }

    //Catch-all case: assume real unless specified otherwise
    template<typename T>
    bool
    operator()(const T& d) { return false; }
    };

bool
isComplex(const IQTensor& T)
    {
    return applyFunc<IsComplex>(T.data());
    }

IQTensor
combiner(std::vector<IQIndex> inds)
    {
    using QNm = pair<QN,long>;

    if(inds.empty()) Error("No indices passed to combiner");
    auto r = inds.size();
    auto type = inds.front().type();
    auto dir = inds.front().dir();

    //create combined index
    detail::GCounter C(r);
    long nsubblocks = 1;
    for(size_t j = 0; j < r; ++j)
        {
        C.setInd(j,1,inds[j].nindex());
        nsubblocks *= inds[j].nindex();
        }
    vector<QNm> qns(nsubblocks);
    for(auto& qm : qns)
        {
        assert(C.notDone());
        qm.second = 1;
        for(size_t j = 0; j < r; ++j)
            {
            qm.first += inds[j].qn(C.i[j])*inds[j].dir();
            qm.second *= inds[j].index(C.i[j]).m();
            }
        ++C;
        }
    //Sort qns by quantum number; there will be duplicates in general:
    auto qnless = [](const QNm& qm1, const QNm& qm2) { return qm1.first < qm2.first; };
    std::sort(qns.begin(),qns.end(),qnless);

    //Count number of sectors
    auto currqn = qns.front().first;
    long sectors = 1;
    for(const auto& qm : qns)
        {
        if(qm.first != currqn)
            {
            ++sectors;
            currqn = qm.first;
            }
        }
    IQIndex::storage indqn(sectors);
    currqn = qns.front().first;
    long currm = 0,
         count = 0;
    for(const auto& qm : qns)
        {
        if(qm.first != currqn)
            {
            indqn[count] = IndexQN(Index(nameint("C",count),currm,type),currqn);
            currqn = qm.first;
            ++count;
            currm = qm.second;
            }
        else
            {
            currm += qm.second;
            }
        }
    //Handle last element of indqn:
    indqn[count] = IndexQN(Index(nameint("C",count),currm,type),currqn);
    //Finally make new combined IQIndex and put at front of inds
    vector<IQIndex> newinds(inds.size()+1);
    auto it = newinds.begin();
    *it = IQIndex("cmb",std::move(indqn),dir);
    for(const auto& I : inds)
        {
        ++it;
        *it = dag(I);
        }

    return IQTensor(QN(),IQIndexSet(std::move(newinds)),make_newdata<ITCombiner>(),{1.0});
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

struct IQNormNoScale : RegisterFunc<IQNormNoScale,Real>
    {
    template<typename T>
    Real
    operator()(const IQTData<T>& d) 
        { 
        Real nrm = 0;
        for(const auto& elt : d.data)
            {
            nrm += std::norm(elt);
            }
        return std::sqrt(nrm);
        }
    };

Real
norm(const IQTensor& T)
    {
#ifdef DEBUG
    if(!T) Error("IQTensor is default initialized");
#endif
    return fabs(T.scale().real0()) *
           applyFunc<IQNormNoScale>(T.data());
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

struct PrintIQT : RegisterFunc<PrintIQT>
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
    void
    operator()(const IQTData<T>& d) const;
    };

template<typename T>
void PrintIQT::
operator()(const IQTData<T>& d) const
    {
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "(omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) 
        {
        s_ << "  ";
        detail::printVal(s_,scalefac*d.data.front());
        return;
        }
        
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
        for(auto i : count(rank))
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
            applyFunc<PrintIQT>(T.data(),s,T.scale(),T.inds());
            }
        }
    else
        {
        if(T.inds()) s << T.inds() << "\n";
        s << "(default constructed)\n";
        }
	s << "\\------------------------------------\n\n";
    return s;
    }


} //namespace itensor
