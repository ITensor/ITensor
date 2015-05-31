//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/iqtdata.h"
#include "itensor/itdata/itdata.h"
#include "itensor/iqindex.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/algs.h"
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"

using std::vector;

namespace itensor {

IQTData::BlOf
make_blof(long b, long o)
    {
    IQTData::BlOf B;
    B.block = b;
    B.offset = o;
    return B;
    }

//function object for calling binaryFind
//on offset vectors below
struct compBlock
    {
    using BlOf = typename IQTData::BlOf;
    bool
    operator()(const BlOf& bo1,
               const BlOf& bo2) const
        { return bo1.block < bo2.block; }
    bool
    operator()(const BlOf& bo, long blk) const        
        { return bo.block < blk; }
    bool
    operator()(long blk, const BlOf& bo) const 
        { return blk < bo.block; }
    };

QN
calcDiv(const IQIndexSet& is, const vector<long>& block_ind)
    {
    QN div;
    for(auto i : count(is.r())) { div += is[i].dir()*is[i].qn(1+block_ind[i]); }
    return div;
    }

void
inverseBlockInd(long block,
                const IQIndexSet& is,
                vector<long>& ind)
    {
    auto r = int(ind.size());
    assert(r == is.r());
    for(int j = 0; j < r-1; ++j)
        {
        ind[j] = block % is[j].nindex();
        block = (block-ind[j])/is[j].nindex();
        }
    ind[r-1] = block;
    }

QN
calcDiv(const IQIndexSet& is, const IQTData& D)
    {
#ifdef DEBUG
    if(D.offsets.empty()) Error("Default constructed IQTData in calcDiv");
#endif
    auto b = D.offsets.front().block;
    QN div;
    auto r = long(is.r());
    for(long j = 0; j < r-1; ++j)
        {
        auto& J = is[j];
        auto Ij = b % J.nindex();
        div += J.dir()*J.qn(1+Ij);
        b = (b-Ij)/J.nindex();
        }
    div += is[r-1].dir()*is[r-1].qn(1+b);
    return div;
    }

IQTData::
IQTData(const IQIndexSet& is, 
        const QN& div)
    {
    auto totalsize = updateOffsets(is,div);
    data.assign(totalsize,0);
    }

long IQTData::
updateOffsets(const IQIndexSet& is,
              const QN& div)
    {
    offsets.clear();

    if(is.r()==0)
        {
        offsets.push_back(make_blof(0,0));
        return 1;
        }

    detail::GCounter C(0,is.r()-1,0);
    for(int j = 0; j < is.r(); ++j) 
        C.setInd(j,0,is[j].nindex()-1);

    long totalsize = 0;
    for(; C.notDone(); ++C)
        {
        QN blockqn;
        for(int j = 0; j < is.r(); ++j)
            {
            auto& J = is[j];
            blockqn += J.qn(1+C.i[j])*J.dir();
            }
        if(blockqn == div)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 totm = 1;   //accumulate area of Indices
            for(int j = 0; j < is.r(); ++j)
                {
                auto& J = is[j];
                auto i_j = C.i[j];
                ind += i_j*indstr;
                indstr *= J.nindex();
                totm *= J[i_j].m();
                }
            offsets.push_back(make_blof(ind,totalsize));
            totalsize += totm;
            }
        }
    return totalsize;
    }

long IQTData::
offsetOf(long blkind) const
    {
    auto blk = detail::binaryFind(offsets,blkind,compBlock());
    if(blk) return blk->offset;
    return -1;
    }

Cplx
doTask(GetElt<IQIndex>& G, const IQTData& d)
    {
    auto* pelt = d.getElt(G.is,G.inds);
    if(pelt) return *pelt;
    return 0;
    }

void
doTask(SetElt<Real,IQIndex>& S, IQTData& d)
    {
    auto* pelt = d.getElt(S.is,S.inds);
    if(pelt) *pelt = S.elt;
    else     Error("Setting IQTensor element non-zero would violate its symmetry.");
    }

void
doTask(MultReal& M, IQTData& d)
    {
    //use BLAS algorithm?
    for(auto& elt : d.data)
        elt *= M.r;
    }


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

    Con.computeScalefac(C.data);
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


void
doTask(Conj, const IQTData& d) { }


Real
doTask(const NormNoScale<IQIndex>& N, const IQTData& d) 
    { 
    Real nrm = 0;
    for(auto& elt : d.data)
        {
        nrm += elt*elt;
        }
    return std::sqrt(nrm);
    }


void
doTask(PrintIT<IQIndex>& P, const IQTData& d)
    {
    P.s << "{Data size = " << d.data.size() << "}\n\n";
    Real scalefac = 1.0;
    if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    else P.s << "(omitting too large scale factor)\n";

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.printVal(scalefac*d.data.front());
        return;
        }
        
    vector<long> block(rank,0);
    auto blockIndex = [&block,&P](long i)->Index { return (P.is[i])[block[i]]; };

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
            if(i > 0) P.s << ", ";
            P.s << blockIndex(i) << "<" << P.is[i].dir() << ">";
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
                P.s << "(";
                for(auto ii = C.i.mini(); ii <= C.i.maxi(); ++ii)
                    {
                    P.s << (1+C.i(ii));
                    if(ii < C.i.maxi()) P.s << ",";
                    }
                P.s << ") ";
                P.printVal(val);
                }
            }
        }

    }

void
doTask(Write& W, const IQTData& d)
    {
    W.writeType(StorageType::IQTData,d); 
    }

} //namespace itensor

