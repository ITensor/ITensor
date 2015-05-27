//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/count.h"
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/detail/printing.h"
#include "itensor/tensor/contract.h"
#include "itensor/iqtensor.h"

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


template<>
IQTensor::
ITensorT(Complex val) 
    :
    scale_(1.)
    { 
    if(val.imag()==0)
        store_ = make_shared<ITDataType<ITDiag<Real>>>(val.real());
    else
        store_ = make_shared<ITDataType<ITDiag<Complex>>>(val);
    }

//IQTensor::
//IQTensor(const QN& q, vector<IQIndex>&& iqinds) 
//	: 
//    is_(move(iqinds)),
//    store_(make_shared<IQTData>(is_,q)),
//    div_(q),
//    scale_(1.)
//	{ }

//IQTensor::
//IQTensor(const QN& q,
//         IQIndexSet&& iset,
//         NewData nd,
//         LogNumber scale)
//    :
//    is_(move(iset)),
//    store_(move(nd)),
//    div_(q),
//    scale_(scale)
//    {
//    }


void
doTask(const PlusEQ<IQIndex>& P,
       IQTData& A,
       const IQTData& B)
    {
#ifdef DEBUG
    if(A.data.size() != B.data.size()) Error("Mismatched sizes in plusEq");
#endif
    if(!P.hasPerm())
        {
        daxpy_wrapper(A.data.size(),P.fac,B.data.data(),1,A.data.data(),1);
        }
    else
        {
        auto r = P.is1().r();
        vector<long> Ablock(r,0),
                     Bblock(r,0);
        Range Arange,
              Brange;
        for(const auto& aio : A.offsets)
            {
            inverseBlockInd(aio.block,P.is1(),Ablock);
            for(int i = 0; i < r; ++i)
                Bblock[i] = Ablock[P.perm().dest(i)];
            Arange.init(make_indexdim(P.is1(),Ablock));
            Brange.init(make_indexdim(P.is2(),Bblock));
            auto* bblock = B.getBlock(P.is2(),Bblock);

            auto aref = makeTensorRef(A.data.data()+aio.offset,Arange);
            auto bref = makeTensorRef(bblock,Brange);
            auto add = [f=P.fac](Real& r1, Real r2) { r1 += f*r2; };
            permute(bref,P.perm(),aref,add);
            }
        }
    }

IQTensor&
operator+=(IQTensor& A, const IQTensor& B)
    {
    if(!A) Error("Calling += on default constructed IQTensor");
    if(!B) Error("Right-hand-side of IQTensor += is default constructed");
    if(&A == &B) return operator*=(A,2.);
    if(A.scale().isZero()) return A.operator=(B);

    Permutation P(A.inds().size());
    try {
        calc_permutation(B.inds(),A.inds(),P);
        }
    catch(const ITError& e)
        {
        Print(A);
        Print(B);
        Error("IQTensor::operator+=: different IQIndex structure");
        }

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
        doTask(PlusEQ<IQIndex>{scalefac},A.store(),B.store());
        }
    else
        {
        doTask(PlusEQ<IQIndex>{P,A.inds(),B.inds(),scalefac},A.store(),B.store());
        }

    return A;
    }

IQTensor&
operator-=(IQTensor& A, const IQTensor& B)
    { 
    if(&A == &B) return operator*=(A,0);
    A.scale().negate();
    operator+=(A,B);
    A.scale().negate();
    return A;
    }

void
permuteIQ(const Permutation& P,
          const IQIndexSet& Ais,
          const IQTData& dA,
          IQIndexSet& Bis,
          IQTData& dB)
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
    dB = std::move(IQTData(Bis,calcDiv(Ais,dA)));

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

        auto* bblock = dB.getBlock(Bis,Bblock);
        auto aref = makeTensorRef(dA.data.data()+aio.offset,Arange);
        auto bref = makeTensorRef(bblock,Brange);
        permute(aref,P,bref);
        }
    }

//class QContract : public RegisterFunc<QContract>
//    {
//    const Label &Aind_,
//                &Bind_;
//
//    const IQIndexSet &Ais_,
//                     &Bis_;
//    QN Cdiv_;
//    IQIndexSet Nis_;
//    Real scalefac_ = -1;

void
doTask(Contract<IQIndex>& Con,
       const IQTData& A,
       const IQTData& B,
       ManagePtr& mp)
    {
    //compute new index set (Con.Nis):
    contractIS(Con.Lis,Con.Lind,Con.Ris,Con.Rind,Con.Nis,true);

    auto Cdiv = calcDiv(Con.Lis,A)+calcDiv(Con.Ris,B);

    //Allocate storage for C
    auto nd = mp.makeNewData<IQTData>(Con.Nis,Cdiv);
    auto& C = *nd;

    auto rA = Con.Lis.r(),
         rB = Con.Ris.r(),
         rC = Con.Nis.r();

    Label AtoB(rA,-1),
          AtoC(rA,-1),
          BtoC(rB,-1);
    Label Cind(rC,0);
    for(auto ic : count(rC))
        {
        auto j = findindex(Con.Lis,Con.Nis[ic]);
        if(j >= 0)
            {
            Cind[ic] = Con.Lind[j];
            AtoC[j] = ic;
            }
        else
            {
            j = findindex(Con.Ris,Con.Nis[ic]);
            Cind[ic] = Con.Rind[j];
            BtoC[j] = ic;
            }
        }
    for(int ia = 0; ia < rA; ++ia)
    for(int ib = 0; ib < rB; ++ib)
        if(Con.Lind[ia] == Con.Rind[ib])
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
        inverseBlockInd(aio.block,Con.Lis,Ablock);
        //Reset couB to run over indices of B (at first)
        couB.reset();
        for(int ib = 0; ib < rB; ++ib)
            couB.setInd(ib,0,Con.Ris[ib].nindex()-1);
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
            auto* bblock = B.getBlock(Con.Ris,couB.i);
            if(!bblock) continue;

            //Finish making Cblock index array
            for(int ib = 0; ib < rB; ++ib)
                if(BtoC[ib] != -1) Cblock[BtoC[ib]] = couB.i[ib];

            auto* cblock = C.getBlock(Con.Nis,Cblock);
            assert(cblock != nullptr);

            //Construct range objects for aref,bref,cref
            //using IndexDim helper objects
            Arange.init(make_indexdim(Con.Lis,Ablock));
            Brange.init(make_indexdim(Con.Ris,couB.i));
            Crange.init(make_indexdim(Con.Nis,Cblock));

            //"Wire up" TensorRef's pointing to blocks of A,B, and C
            //we are working with
            auto aref = makeTensorRef(A.data.data()+aio.offset,Arange),
                 bref = makeTensorRef(bblock,Brange);
            auto cref= makeTensorRef(cblock,Crange);

            //Compute aref*bref=cref
            contract(aref,Con.Lind,bref,Con.Rind,cref,Cind);

            } //for couB
        } //for A.offsets

    //Compute new scalefac_ from C.data
    Con.scalefac = 0;
    for(auto elt : C.data) Con.scalefac += elt*elt;
    Con.scalefac = std::sqrt(Con.scalefac);
    //Rescale C by scalefac
    if(Con.scalefac != 0)
        {
        for(auto& elt : C.data) elt /= Con.scalefac;
        }
    }

void
combine(const IQTData& d,
        const IQIndexSet& dis,
        const IQIndexSet& Cis,
        IQIndexSet& Nis,
        ManagePtr& mp,
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

    IQTData nd;
    if(P) 
        {
        permuteIQ(P,dis,d,Nis,nd);
        }
    auto& Pis = (P ? Nis : dis);

    if(jc == 0) //has cind at front, we are "uncombining"
        {
        auto newr = Pis.r()+Cis.r()-2;
        auto offset = Cis.r()-1;
        vector<IQIndex> newind(newr);
        for(auto j : count(offset))
            newind.at(j) = Cis[1+j];
        for(auto j : count(Pis.r()-1))
            newind.at(offset+j) = Pis[1+j];
        Nis = IQIndexSet(move(newind));
        }
    else //we are "combining"
        {
        auto newr = Pis.r()-Cis.r()+2;
        auto offset = Cis.r()-2;
        vector<IQIndex> newind(newr);
        newind.front() = cind;
        for(auto j : count(1,newr))
            newind.at(j) = Pis[offset+j];
        Nis = IQIndexSet(move(newind));
        }

    //Only need to modify d if Cis.r() > 2.
    //If Cis.r()==2 just swapping one index for another
    if(Cis.r() > 2)
        {
        IQTData* p = nullptr;
        if(nd)
            {
            p = &nd;
            }
        else if(own_data) 
            {
            p = mp.modifyData(d);
            }
        else
            {
            nd = d;
            p = &nd;
            }
        auto div = calcDiv(dis,d);
        p->updateOffsets(Nis,div);
        }

    if(nd) mp.makeNewData<IQTData>(std::move(nd));
    }

void
doTask(Contract<IQIndex>& C,
       const IQTData& d,
       const ITCombiner& cmb,
       ManagePtr& mp)
    {
    combine(d,C.Lis,C.Ris,C.Nis,mp,true);
    }

void
doTask(Contract<IQIndex>& C,
       const ITCombiner& cmb,
       const IQTData& d,
       ManagePtr& mp)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,mp,false);
    }


IQTensor&
operator*=(IQTensor& A, const IQTensor& B)
    {
    if(!A || !B)
        Error("Default constructed IQTensor in product");

    if(&A == &B)
        {
        A = IQTensor(sqr(norm(A)));
        return A;
        }

    auto& Lis = A.inds();
    auto& Ris = B.inds();

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

    auto nstore = A.store();
    auto C = doTask(Contract<IQIndex>{Lis,Lind,Ris,Rind},nstore,B.store());

    auto nscale = A.scale()*B.scale();
    if(!std::isnan(C.scalefac)) nscale *= C.scalefac;

#ifdef DEBUG
    //Check for duplicate indices
    detail::check(C.Nis);
#endif

    A = IQTensor(C.Nis,std::move(nstore),nscale);

    return A;
    }

struct AddITensor
    {
    const QN& tdiv;
    const IQIndexSet& iqis;
    const IndexSet& is;
    const vector<long>& block_ind;
    const Permutation& P;
    Real fac = 0;
    AddITensor(const QN& tdiv_,
               const IQIndexSet& iqis_,
               const IndexSet& is_,
               const vector<long>& block_ind_,
               const Permutation& P_,
               Real scalefac_)
        :
        tdiv(tdiv_),
        iqis(iqis_),
        is(is_),
        block_ind(block_ind_),
        P(P_),
        fac(scalefac_)
        { }
    };

void
doTask(AddITensor& A, IQTData& d, const ITReal& t)
    {
    auto ddiv = calcDiv(A.iqis,d);
    if(ddiv != A.tdiv) Error("IQTensor+=ITensor, ITensor has incompatible QN flux/divergence");
    Range drange;
    drange.init(make_indexdim(A.iqis,A.block_ind));
    auto* dblock = d.getBlock(A.iqis,A.block_ind);

    auto dref = makeTensorRef(dblock,drange);
    auto tref = makeTensorRef(t.data(),A.is);
    auto add = [f=A.fac](Real& r1, Real r2) { r1 += f*r2; };
    permute(tref,A.P,dref,add);
    }

IQTensor&
operator+=(IQTensor& T, const ITensor& t)
    {
    if(!t) Error("IQTensor+=ITensor: r.h.s. ITensor is default constructed");
    if(!T.inds()) Error("Calling IQTensor+= on default constructed ITensor");
    auto rank = T.r();
#ifdef DEBUG
    if(t.r() != rank) Error("Mismatched number of indices in IQTensor+=ITensor");
#endif

    Permutation P(rank);
    vector<long> block_ind(rank);
    for(auto i : count(rank))
    for(auto I : count(rank))
        {
        auto j = findindex(T.inds()[I],t.inds()[i]);
        if(j > 0)
            {
            block_ind[I] = (j-1);
            P.setFromTo(i,I);
            break;
            }
        }

    auto tdiv = calcDiv(T.inds(),block_ind);

    if(!T.store()) 
        {
        //allocate data to add this ITensor into
        if(!isComplex(t)) T.store() = make_shared<ITDataType<IQTData>>(T.inds(),tdiv);
        else              Error("Initializing complex IQTensor in +=ITensor not yet implemented");
        }

    Real scalefac = 1;
    if(T.scale().magnitudeLessThan(t.scale())) T.scaleTo(t.scale()); 
    else                                       scalefac = (t.scale()/T.scale()).real();

    doTask(AddITensor{tdiv,T.inds(),t.inds(),block_ind,P,scalefac},T.store(),t.store());

    return T;
    }

template<>
IQTensor& IQTensor::
conj()
    {
    doTask(Conj{},store_);
    return *this;
    }

template<>
IQTensor& IQTensor::
dag()
    {
    is_.dag();
    doTask(Conj{},store_);
    return *this;
    }

void
doTask(MultReal& M, IQTData& d)
    {
    //use BLAS algorithm?
    for(auto& elt : d.data)
        elt *= M.r;
    }

template<>
void IQTensor::
scaleTo(const LogNumber& newscale)
    {
    if(scale_ == newscale) return;
    if(newscale.sign() == 0) Error("Trying to scale an IQTensor to a 0 scale");
    scale_ /= newscale;
    doTask(MultReal{scale_.real0()},store_);
    scale_ = newscale;
    }

struct ToITensor
    {
    const IQIndexSet& is;
    const LogNumber& scale;

    ToITensor(const IQIndexSet& is_,
              const LogNumber& scale_)
        :
        is(is_),
        scale(scale_)
        { }
    };

ITensor
doTask(ToITensor& T, const IQTData& d)
    {
    auto r = T.is.r();
    auto nd = ITReal(area(T.is),0);
    auto *pd = d.data.data();
    auto *pn = nd.data();
    vector<long> block(r,0);
    detail::GCounter C(r);
    for(const auto& io : d.offsets)
        {
        inverseBlockInd(io.block,T.is,block);
        for(long j = 0; j < r; ++j)
            {
            long start = 0;
            for(long b = 0; b < block[j]; ++b)
                start += T.is[j][b].m();
            C.setInd(j,start,start+T.is[j][block[j]].m()-1);
            }
        //TODO: need to make a Range/TensoRef iterator
        //to rewrite the following code more efficiently
        for(; C.notDone(); ++C)
            {
            pn[ind(T.is,C.i)] = pd[io.offset+C.ind];
            }
        }
    vector<Index> inds(r);
    for(long j = 0; j < r; ++j) inds[j] = T.is[j];
    return ITensor(IndexSet(std::move(inds)),std::move(nd),T.scale);
    }

ITensor
toITensor(const IQTensor& T)
    {
    return doTask<ITensor>(ToITensor{T.inds(),T.scale()},T.store());
    }

bool
doTask(IsComplex,const IQTData& d) { return false; }

bool
isComplex(const IQTensor& T)
    {
    return doTask<bool>(IsComplex{},T.store());
    }

struct CalcDiv 
    { 
    const IQIndexSet& is;
    CalcDiv(const IQIndexSet& is_) : is(is_) { }
    };

QN
doTask(const CalcDiv& C,const IQTData& d)
    {
    return calcDiv(C.is,d);
    }

QN
div(const IQTensor& T) 
    { 
    if(!T) Error("div(IQTensor) not defined for unallocated IQTensor");
    return doTask<QN>(CalcDiv{T.inds()},T.store());
    }

IQTensor
combiner(std::vector<IQIndex> inds,
         const Args& args)
    {
    using QNm = pair<QN,long>;
    if(inds.empty()) Error("No indices passed to combiner");
    auto cname = args.getString("IndexName","cmb");
    auto itype = getIndexType(args,"IndexType",inds.front().type());
    auto r = inds.size();
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
            indqn[count] = IndexQN(Index(nameint("C",count),currm,itype),currqn);
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
    indqn[count] = IndexQN(Index(nameint("C",count),currm,itype),currqn);
    //Finally make new combined IQIndex and put at front of inds
    vector<IQIndex> newinds(inds.size()+1);
    auto it = newinds.begin();
    *it = IQIndex(cname,std::move(indqn),dir);
    for(const auto& I : inds)
        {
        ++it;
        *it = dag(I);
        }

    return IQTensor(IQIndexSet{std::move(newinds)},ITCombiner());
    }

IQIndex
findIQInd(const IQTensor& T, const Index& i)
    {
    for(const IQIndex& J : T.inds())
        if(hasindex(J,i)) return J;
    Print(T.inds());
    Print(i);
    throw ITError("Index i not found in any of T's IQIndices");
    return IQIndex{};
    }


Arrow
dir(const IQTensor& T, const IQIndex& I)
	{
    for(const IQIndex& J : T.inds())
        if(I == J) return J.dir();
    Error("dir: IQIndex not found");
    return Out;
	}

Real
doTask(const NormNoScale<IQIndex>& N, const IQTData& d) 
    { 
    Real nrm = 0;
    for(const auto& elt : d.data)
        {
        nrm += std::norm(elt);
        }
    return std::sqrt(nrm);
    }

Real
norm(const IQTensor& T)
    {
#ifdef DEBUG
    if(!T) Error("IQTensor is default initialized");
#endif
    return fabs(T.scale().real0()) *
           doTask<Real>(NormNoScale<IQIndex>{T.inds()},T.store());
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


void
doTask(PrintIT<IQIndex>& P, const IQTData& d)
    {
    Real scalefac = 1.0;
    if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    else P.s << "(omitting too large scale factor)\n";

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        detail::printVal(P.s,scalefac*d.data.front());
        return;
        }
        
    vector<long> block(rank,0);
    auto blockIndex = [&block,&P](long i)->const Index& { return (P.is[i])[block[i]]; };

    Range brange;
    detail::GCounter C(rank);
    for(const auto& io : d.offsets)
        {
        //Determine block indices (where in the IQIndex space
        //this non-zero block is located)
        inverseBlockInd(io.block,P.is,block);
        //Print Indices of this block
        for(auto i : count(rank))
            {
            P.s << blockIndex(i) << "<" << P.is[i].dir() << "> ";
            }
        P.s << "\n";
        //Wire up GCounter with appropriate dims
        C.reset();
        for(int i = 0; i < rank; ++i)
            C.setInd(i,0,blockIndex(i).m()-1);
        for(auto os = io.offset; C.notDone(); ++C, ++os)
            {
            auto val = scalefac*d.data[os];
            if(std::norm(val) > Global::printScale())
                {
                P.s << "  (";
                for(auto ii = C.i.mini(); ii <= C.i.maxi(); ++ii)
                    {
                    P.s << (1+C.i(ii));
                    if(ii < C.i.maxi()) P.s << ",";
                    }
                P.s << ") ";
                detail::printVal(P.s,val);
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
        bool ff_set = (std::ios::floatfield & s.flags()) != 0;

        if(ff_set || Global::printdat())
            {
            s << "\n";
            doTask(PrintIT<IQIndex>{s,T.scale(),T.inds(),true},T.store());
            }
        }
    else
        {
        if(T.inds()) s << T.inds() << "\n";
        s << "(storage not allocated)\n";
        }
	s << "\\------------------------------------\n\n";
    return s;
    }


} //namespace itensor
