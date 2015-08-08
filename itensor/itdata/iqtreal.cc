//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/iqtreal.h"
#include "itensor/itdata/itdata.h"
#include "itensor/iqindex.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/algs.h"
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"

using std::vector;
using std::move;

namespace itensor {


//function object for calling binaryFind
//on offset vectors below
struct compBlock
    {
    using BlOf = typename IQTReal::BlOf;
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
calcDiv(IQIndexSet const& is, Label const& block_ind)
    {
    QN div;
    for(auto i : count(is.r())) { div += is[i].dir()*is[i].qn(1+block_ind[i]); }
    return div;
    }

QN
doTask(CalcDiv const& C,
       IQTReal const& D)
    {
#ifdef DEBUG
    if(D.offsets.empty()) Error("Default constructed IQTReal in doTask(CalcDiv,IQTReal)");
#endif
    auto b = D.offsets.front().block;
    Label block_ind(C.is.r());
    inverseBlockInd(b,C.is,block_ind);
    return calcDiv(C.is,block_ind);

    //auto r = long(is.r());
    //if(r==0) return div;

    //for(long j = 0; j < r-1; ++j)
    //    {
    //    auto& J = is[j];
    //    auto Ij = b % J.nindex();
    //    div += J.dir()*J.qn(1+Ij);
    //    b = (b-Ij)/J.nindex();
    //    }
    //div += is[r-1].dir()*is[r-1].qn(1+b);
    //return div;
    }

IQTReal::
IQTReal(const IQIndexSet& is, 
        const QN& div)
    {
    auto totalsize = updateOffsets(is,div);
    store.assign(totalsize,0);
    }

long IQTReal::
updateOffsets(const IQIndexSet& is,
              const QN& div)
    {
    offsets.clear();

    if(is.r()==0)
        {
        offsets.push_back(make_blof(0,0));
        return 1;
        }

    //Set up counter over all blocks
    auto C = detail::GCounter(is.r());
    for(auto j : count(is.r()))
        C.setRange(j,0,is[j].nindex()-1);

    long totalsize = 0;
    for(; C.notDone(); ++C)
        {
        QN blockqn;
        for(auto j : count(is.r()))
            {
            auto& J = is[j];
            blockqn += J.qn(1+C[j])*J.dir();
            }
        if(blockqn == div)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 totm = 1;   //accumulate area of Indices
            for(auto j : count(is.r()))
                {
                auto& J = is[j];
                auto i_j = C[j];
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

long
offsetOf(std::vector<IQTReal::BlOf> const& offsets,
         long blockind)
    {
    auto blk = detail::binaryFind(offsets,blockind,compBlock());
    if(blk) return blk->offset;
    return -1;
    }

Cplx
doTask(GetElt<IQIndex>& G, const IQTReal& d)
    {
    auto* pelt = d.getElt(G.is,G.inds);
    if(pelt) return *pelt;
    return 0;
    }

void
doTask(SetElt<Real,IQIndex>& S, IQTReal& d)
    {
    auto* pelt = d.getElt(S.is,S.inds);
    if(pelt) *pelt = S.elt;
    else     Error("Setting IQTensor element non-zero would violate its symmetry.");
    }

void
doTask(MultReal& M, IQTReal& d)
    {
    //use BLAS algorithm?
    for(auto& elt : d.store)
        elt *= M.r;
    }


void
doTask(PlusEQ<IQIndex> const& P,
       IQTReal & A,
       IQTReal const& B)
    {
#ifdef DEBUG
    if(A.store.size() != B.store.size()) Error("Mismatched sizes in plusEq");
#endif
    if(!P.hasPerm())
        {
        daxpy_wrapper(A.store.size(),P.fac,B.data(),1,A.data(),1);
        }
    else
        {
        auto r = P.is1().r();
        Label Ablock(r,0),
              Bblock(r,0);
        Range Arange,
              Brange;
        for(auto& aio : A.offsets)
            {
            inverseBlockInd(aio.block,P.is1(),Ablock);
            for(int i = 0; i < r; ++i)
                Bblock[i] = Ablock[P.perm().dest(i)];
            Arange.init(make_indexdim(P.is1(),Ablock));
            Brange.init(make_indexdim(P.is2(),Bblock));
            auto* bblock = getBlock(B,P.is2(),Bblock);

            auto aref = makeTenRef(A.data()+aio.offset,Arange);
            auto bref = makeTenRef(bblock,Brange);
            auto add = [f=P.fac](Real& r1, Real r2) { r1 += f*r2; };
            permute(bref,P.perm(),aref,add);
            }
        }
    }


void
doTask(Contract<IQIndex>& Con,
       IQTReal const& A,
       IQTReal const& B,
       ManageStore& m)
    {
    Label Lind,
          Rind;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),Lind,Rind);
    //compute new index set (Con.Nis):
    Label Cind;
    contractIS(Con.Lis,Lind,Con.Ris,Rind,Con.Nis,Cind,true);

    auto Cdiv = doTask(CalcDiv{Con.Lis},A)+doTask(CalcDiv{Con.Ris},B);

    //Allocate storage for C
    auto nd = m.makeNewData<IQTReal>(Con.Nis,Cdiv);
    auto& C = *nd;

    //Function to execute for each pair of
    //contracted blocks of A and B
    auto do_contract = 
        [&Con,&Lind,&Rind,&Cind]
        (const Real *ablock, Label const& Ablockind,
         const Real *bblock, Label const& Bblockind,
               Real *cblock, Label const& Cblockind)
        {
        Range Arange,
              Brange,
              Crange;
        //Construct range objects for aref,bref,cref
        //using IndexDim helper objects
        Arange.init(make_indexdim(Con.Lis,Ablockind));
        Brange.init(make_indexdim(Con.Ris,Bblockind));
        Crange.init(make_indexdim(Con.Nis,Cblockind));

        //"Wire up" TensorRef's pointing to blocks of A,B, and C
        //we are working with
        auto aref = makeTenRef(ablock,Arange),
             bref = makeTenRef(bblock,Brange);
        auto cref = makeTenRef(cblock,Crange);

        //Compute cref=aref*bref
        contract(aref,Lind,bref,Rind,cref,Cind);
        };

    loopContractedBlocks(A,Con.Lis,
                         B,Con.Ris,
                         C,Con.Nis,
                         do_contract);

    Con.computeScalefac(C.store);
    }

void
permuteIQ(const Permutation& P,
          const IQIndexSet& Ais,
          const IQTReal& dA,
          IQIndexSet& Bis,
          IQTReal& dB)
    {
#ifdef DEBUG
    if(isTrivial(P)) Error("Calling permuteIQ for trivial Permutation");
#endif
    auto r = Ais.r();
    auto bind = IQIndexSetBuilder(r);
    for(auto i : count(r))
        {
        bind.setExtent(P.dest(i),Ais[i]);
        }
    Bis = IQIndexSet(bind);
    dB = IQTReal(Bis,doTask(CalcDiv{Ais},dA));

    Label Ablock(r,-1),
          Bblock(r,-1);
    Range Arange,
          Brange;
    if(Global::debug1()) println("P = ",P);
    for(auto aio : dA.offsets)
        {
        //Compute bi, new block index of blk
        inverseBlockInd(aio.block,Ais,Ablock);
        for(auto j : index(Ablock)) 
            Bblock.at(P.dest(j)) = Ablock[j];
        Arange.init(make_indexdim(Ais,Ablock));
        Brange.init(make_indexdim(Bis,Bblock));

        auto* bblock = getBlock(dB,Bis,Bblock);
        auto aref = makeTenRef(dA.data()+aio.offset,Arange);
        auto bref = makeTenRef(bblock,Brange);
        permute(aref,P,bref);
        }
    }

IQIndexSet
replaceInd(const IQIndexSet& is,
           long loc,
           const IQIndex& replacement)
    {
    auto newind = IQIndexSetBuilder(is.r());
    long i = 0;
    for(long j = 0; j < loc; ++j)
        newind.setExtent(i++,is[j]);
    newind.setExtent(i++,replacement);
    for(decltype(is.r()) j = loc+1; j < is.r(); ++j)
        newind.setExtent(i++,is[j]);
    return IQIndexSet(newind);
    }

void
condense(IQIndexSet const& Cis,
         IQIndexSet const& dis,
         IQTReal & d)
    {
    auto r = dis.r();
    auto is_cmb = InfArray<int,10ul>(r,0);
    for(auto i : count(r))
        if(hasindex(Cis,dis[i])) is_cmb[i] = 1;

    //Loop over non-zero blocks
    auto block = Label(r);
    for(auto io : d.offsets)
        {
        inverseBlockInd(io.block,dis,block);
        }
    }

void
combine(IQTReal const& d,
        IQIndexSet const& dis,
        IQIndexSet const& Cis,
        IQIndexSet & Nis,
        ManageStore & m,
        bool own_data)
    {
    //cind is special "combined index"
    auto const& cind = Cis[0];
    //check if d has combined index i.e. we are "uncombining"
    auto jc = findindex(dis,cind);
    auto combining = (jc < 0);
    auto uncombining = not combining;

    IQTReal* pd = nullptr;
    //Call this if necessary to modify the data
    auto setPtrData = [&]()
        {
        if(pd) return;
        if(own_data) pd = m.modifyData(d);
        else         pd = m.makeNewData<IQTReal>(d);
        };

    Permutation P;

    if(Cis.r() == 2) //treat rank 2 combiner specially
        {
        if(uncombining)
            {
            //Has cind, replace with Cis[1]
            Nis = replaceInd(dis,jc,Cis[1]);
            }
        else //combining
            {
            //Has Cis[1], replace with cind
            auto ju = findindex(dis,Cis[1]);
            if(ju < 0)
                {
                println("IQIndexSet of regular IQTensor =\n",dis);
                println("IQIndexSet of combiner/delta =\n",Cis);
                println("Missing IQIndex: ",Cis[1]);
                Error("IQCombiner: missing IQIndex");
                }
            Nis = replaceInd(dis,ju,cind);
            }
        if(!own_data) m.assignPointerRtoL();
        return;
        }
    else if(uncombining && jc != 0) //we are uncombining, but cind not at front
        {
        P = Permutation(dis.r());
        P.setFromTo(jc,0);
        long ni = 1;
        for(auto j : count(dis.r()))
            if(j != jc) P.setFromTo(j,ni++);
        jc = 0;
        }
    else if(combining) //we are combining, set up Permutation P
        {
        //check locations of Cis[1], Cis[2], ...
        //Check if Cis[1],Cis[2],... are grouped together (contiguous);
        //all at front; and in same order as on combiner
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
                    println("IQIndexSet of regular IQTensor =\n",dis);
                    println("IQIndexSet of combiner/delta =\n",Cis);
                    println("Missing IQIndex: ",Cis[c]);
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

    if(P) 
        {
        setPtrData(); //sets pd
        permuteIQ(P,dis,d,Nis,*pd);
        }
    //Pis means 'permuted' index set
    auto& Pis = (P ? Nis : dis);

    if(uncombining)
        {
        auto newr = Pis.r()+Cis.r()-2;
        auto offset = Cis.r()-1;
        auto newind = IQIndexSetBuilder(newr);
        for(auto j : count(offset))
            newind.setExtent(j,Cis[1+j]);
        for(auto j : count(Pis.r()-1))
            newind.setExtent(offset+j,Pis[1+j]);
        Nis = IQIndexSet(newind);
        }
    else //combining
        {
        auto newr = Pis.r()-Cis.r()+2;
        auto newind = IQIndexSetBuilder(newr);
        newind.setExtent(0,cind);
        for(auto j : count(1,newr)) newind.setExtent(j,Pis[Cis.r()-2+j]);
        setPtrData();
        condense(Cis,Pis,*pd);
        Nis = IQIndexSet(newind);
        }

    //Only need to modify d if Cis.r() > 2.
    //If Cis.r()==2 just swapping one index for another
    if(Cis.r() > 2)
        {
        setPtrData();
        auto div = doTask(CalcDiv{dis},d);
        pd->updateOffsets(Nis,div);
        }
    }

void
doTask(Contract<IQIndex>& C,
       const IQTReal& d,
       const ITCombiner& cmb,
       ManageStore& m)
    {
    combine(d,C.Lis,C.Ris,C.Nis,m,true);
    }

void
doTask(Contract<IQIndex>& C,
       const ITCombiner& cmb,
       const IQTReal& d,
       ManageStore& m)
    { 
    combine(d,C.Ris,C.Lis,C.Nis,m,false);
    }


void
doTask(Conj, const IQTReal& d) { }


Real
doTask(NormNoScale, const IQTReal& d) 
    { 
    Real nrm = 0;
    for(auto& elt : d.store)
        {
        nrm += elt*elt;
        }
    return std::sqrt(nrm);
    }


void
doTask(PrintIT<IQIndex>& P, const IQTReal& d)
    {
    P.s << "IQTReal {" << d.offsets.size() << " blocks; Data size = " << d.store.size() << "}\n\n";
    Real scalefac = 1.0;
    if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    else P.s << "(omitting too large scale factor)\n";

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.printVal(scalefac*d.store.front());
        return;
        }
        
    Label block(rank,0);
    auto blockIndex = [&block,&P](long i)->Index { return (P.is[i])[block[i]]; };

    Range brange;
    auto C = detail::GCounter(rank);
    for(const auto& io : d.offsets)
        {
        bool indices_printed = false;
        //Determine block indices (where in the IQIndex space
        //this non-zero block is located)
        inverseBlockInd(io.block,P.is,block);

        Label boff(rank,0);
        for(auto i : count(rank))
            {
            for(auto j : count(block[i]))
                boff[i] += P.is[i][j].m();
            }

        //Wire up GCounter with appropriate dims
        C.reset();
        for(decltype(rank) i = 0; i < rank; ++i)
            C.setRange(i,0,blockIndex(i).m()-1);
        for(auto os = io.offset; C.notDone(); ++C, ++os)
            {
            auto val = scalefac*d.store[os];
            if(std::norm(val) > Global::printScale())
                {
                if(!indices_printed)
                    {
                    indices_printed = true;
                    //Print Indices of this block
                    for(auto i : count(rank))
                        {
                        if(i > 0) P.s << ", ";
                        P.s << blockIndex(i) << "<" << P.is[i].dir() << ">";
                        }
                    P.s << "\n";
                    }
                P.s << "(";
                for(auto ii : count(rank))
                    {
                    P.s << (1+boff[ii]+C[ii]);
                    if(1+ii != rank) P.s << ",";
                    }
                P.s << ") ";

                //P.s << "[";
                //for(auto ii : count(rank))
                //    {
                //    P.s << (1+C[ii]);
                //    if(1+ii != rank) P.s << ",";
                //    }
                //P.s << "] ";

                P.printVal(val);
                }
            }
        }
    }


void
doTask(Write& W, const IQTReal& d)
    {
    W.writeType(StorageType::IQTReal,d); 
    }

} //namespace itensor

