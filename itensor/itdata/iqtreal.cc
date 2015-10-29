//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/count.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/algs.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"
#include "itensor/itdata/iqtreal.h"

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
calcDiv(IQIndexSet const& is, 
        Label const& block_ind)
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
    computeBlockInd(b,C.is,block_ind);
    return calcDiv(C.is,block_ind);
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

    //Set up a Range to iterate over all blocks
    auto RB = RangeBuilder(is.r());
    for(auto j : count(is.r()))
        RB.nextIndex(is[j].nindex());

    long totalsize = 0;
    for(auto I : RB.build())
        {
        auto blockqn = QN{};
        for(auto j : count(is.r()))
            {
            auto& J = is[j];
            blockqn += J.qn(1+I[j])*J.dir();
            }
        if(blockqn == div)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 totm = 1;   //accumulate area of Indices
            for(auto j : count(is.r()))
                {
                auto& J = is[j];
                auto i_j = I[j];
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

Cplx
doTask(SumEls<IQIndex>, IQTReal const& d)
    {
    Real s = 0.;
    for(auto& el : d.store) s += el;
    return Cplx(s,0.);
    }

void
doTask(MultReal const& M, IQTReal & d)
    {
    dscal_wrapper(d.store.size(),M.r,d.store.data());
    }


void
doTask(PlusEQ<IQIndex> const& P,
       IQTReal & A,
       IQTReal const& B)
    {
#ifdef DEBUG
    if(A.store.size() != B.store.size()) Error("Mismatched sizes in plusEq");
#endif
    if(isTrivial(P.perm()))
        {
        daxpy_wrapper(A.store.size(),P.fac(),B.data(),1,A.data(),1);
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
            computeBlockInd(aio.block,P.is1(),Ablock);
            for(int i = 0; i < r; ++i)
                Bblock[i] = Ablock[P.perm().dest(i)];
            Arange.init(make_indexdim(P.is1(),Ablock));
            Brange.init(make_indexdim(P.is2(),Bblock));

            auto aref = makeTenRef(A.data(),aio.offset,A.size(),&Arange);

            auto bblock = getBlock(B,P.is2(),Bblock);
            auto bref = TensorRefc(bblock,&Brange);

            //aref += permute(bref,P.perm());
            auto f = P.fac();
            auto add = [f](Real r2, Real& r1) { r1 += f*r2; };
            transform(permute(bref,P.perm()),aref,add);
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
    const bool sortResult = false;
    contractIS(Con.Lis,Lind,Con.Ris,Rind,Con.Nis,Cind,sortResult);

    auto Cdiv = doTask(CalcDiv{Con.Lis},A)+doTask(CalcDiv{Con.Ris},B);

    //Allocate storage for C
    START_TIMER(33)
    auto nd = m.makeNewData<IQTReal>(Con.Nis,Cdiv);
    STOP_TIMER(33)
    auto& C = *nd;

    //Function to execute for each pair of
    //contracted blocks of A and B
    auto do_contract = 
        [&Con,&Lind,&Rind,&Cind]
        (Datac ablock, Label const& Ablockind,
         Datac bblock, Label const& Bblockind,
         Data  cblock, Label const& Cblockind)
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
        auto aref = TensorRefc(ablock,&Arange),
             bref = TensorRefc(bblock,&Brange);
        auto cref = TensorRef(cblock,&Crange);

        //Compute cref += aref*bref
        START_TIMER(2)
        contract(aref,Lind,bref,Rind,cref,Cind,1.,1.);
        STOP_TIMER(2)
        };

    START_TIMER(20)
    loopContractedBlocks(A,Con.Lis,
                         B,Con.Ris,
                         C,Con.Nis,
                         do_contract);
    STOP_TIMER(20)

    START_TIMER(21)
    Con.scalefac = computeScalefac(C);
    STOP_TIMER(21)
    }

void
doTask(NCProd<IQIndex>& P,
       IQTReal const& A,
       IQTReal const& B,
       ManageStore& m)
    {
    auto& Ais = P.Lis;
    auto& Bis = P.Ris;
    auto& Cis = P.Nis;
    auto rA = rank(Ais);
    auto rB = rank(Bis);
    Label Aind,
          Bind,
          Cind;
    computeLabels(Ais,rA,Bis,rB,Aind,Bind);
    ncprod(Ais,Aind,Bis,Bind,Cis,Cind);

    Label BtoA(rA,-1);
    for(auto ia : count(rA))
    for(auto ib : count(rB))
        if(Bis[ib] == Ais[ia])
            {
            BtoA[ib] = ia;
            break;
            }

    auto Cdiv = QN{};
        {
        Cdiv = doTask(CalcDiv{Ais},A);
        auto Ablock_ind = Label(rA);
        computeBlockInd(A.offsets.front().block,Ais,Ablock_ind);
        auto Bblock_ind = Label(rB);
        for(auto& bo : B.offsets)
            {
            computeBlockInd(bo.block,Bis,Bblock_ind);
            bool matchesA = true;
            for(auto n : count(rB))
                {
                if(Bind[n] < 0 && Ablock_ind[BtoA[n]] != Bind[n])
                    {
                    matchesA = false;
                    break;
                    }
                }
            if(matchesA) break;
            }
        //Only account for unique indices of B
        for(auto n : count(rB))
            if(Bind[n] > 0) //unique
                {
                Cdiv += Bis[n].dir()*Bis[n].qn(1+Bblock_ind[n]);
                }
        }

    //Allocate storage for C
    auto& C = *m.makeNewData<IQTReal>(Cis,Cdiv);

    auto do_ncprod = 
        [&P,&Aind,&Bind,&Cind]
        (Datac ablock, Label const& Ablockind,
         Datac bblock, Label const& Bblockind,
         Data  cblock, Label const& Cblockind)
        {
        Range Arange,
              Brange,
              Crange;
        //Construct range objects for aref,bref,cref
        //using IndexDim helper objects
        Arange.init(make_indexdim(P.Lis,Ablockind));
        Brange.init(make_indexdim(P.Ris,Bblockind));
        Crange.init(make_indexdim(P.Nis,Cblockind));

        //"Wire up" TensorRef's pointing to blocks of A,B, and C
        //we are working with
        auto aref = TensorRefc(ablock,&Arange),
             bref = TensorRefc(bblock,&Brange);
        auto cref = TensorRef(cblock,&Crange);

        //Compute cref += aref*bref
        ncprod(aref,Aind,bref,Bind,cref,Cind);
        };

    loopContractedBlocks(A,Ais,
                         B,Bis,
                         C,Cis,
                         do_ncprod);

    P.scalefac = computeScalefac(C);
    }

void
doTask(Conj, IQTReal const& d) { }


Real
doTask(NormNoScale, IQTReal const& d) 
    { 
    return dnrm2_wrapper(d.size(),d.data());
    }


void
doTask(PrintIT<IQIndex>& P, IQTReal const& d)
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
        computeBlockInd(io.block,P.is,block);

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
                        if(i > 0) P.s << " ";
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
doTask(Write& W, IQTReal const& d)
    {
    W.writeType(StorageType::IQTReal,d); 
    }

} //namespace itensor

