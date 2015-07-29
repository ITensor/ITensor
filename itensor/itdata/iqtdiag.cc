//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/iqtdiag.h"
#include "itensor/detail/gcounter.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"

using std::vector;
using std::move;

namespace itensor {

IQTDiag::
IQTDiag(IQIndexSet const& is, 
        QN const& div)
    {
    auto totalsize = updateOffsets(is,div);
    store.assign(totalsize,0);
    }

long IQTDiag::
updateOffsets(IQIndexSet const& is,
              QN const& div)
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
                 minm = std::numeric_limits<long>::max();
            for(int j = 0; j < is.r(); ++j)
                {
                auto& J = is[j];
                auto i_j = C.i[j];
                ind += i_j*indstr;
                indstr *= J.nindex();
                minm = std::min(minm,J[i_j].m());
                }
            offsets.push_back(make_blof(ind,totalsize));
            totalsize += minm;
            }
        }
    return totalsize;
    }

QN
doTask(CalcDiv const& C, IQTDiag const& d)
    {
#ifdef DEBUG
    if(d.offsets.empty()) Error("Default constructed IQTReal in doTask(CalcDiv,IQTReal)");
#endif
    auto b = d.offsets.front().block;
    Label block_ind(C.is.r());
    inverseBlockInd(b,C.is,block_ind);
    return calcDiv(C.is,block_ind);
    }

void
doTask(MultReal & M, IQTDiag & d)
    {
    //use BLAS algorithm?
    for(auto& elt : d.store)
        elt *= M.r;
    }

void
blockDiagDense(IQTDiag const& D,
               IQIndexSet const& Dis,
               Label const& Dind,
               IQTReal const& T,
               IQIndexSet const& Tis,
               Label const& Tind,
               IQIndexSet const& Cis,
               Label const& Cind,
               ManageStore & m)
    {
#ifdef DEBUG
    if(Dis.r() == 0) Error("IQTDiag rank 0 case not handled");
#endif

    bool T_has_uncontracted = false;
    for(auto j : index(Tind)) 
        if(Tind[j] >= 0)
            {
            T_has_uncontracted = true;
            break;
            }

    auto Cdiv = doTask(CalcDiv{Dis},D)+doTask(CalcDiv{Tis},T);

    if(T_has_uncontracted)
        {
        auto *nd = m.makeNewData<IQTReal>(Cis,Cdiv);
        auto& C = *nd;

        auto do_contract =
            [&Dis,&Tis,&Cis,&Dind,&Tind,&Cind]
            (const Real *dblock, Label const& Dblockind,
             const Real *tblock, Label const& Tblockind,
                   Real *cblock, Label const& Cblockind)
            {
            Range Trange,
                  Crange;
            Trange.init(make_indexdim(Tis,Tblockind));
            auto Tref = makeTensorRef(tblock,Trange);
            Crange.init(make_indexdim(Cis,Cblockind));
            auto Cref = makeTensorRef(cblock,Crange);

            auto Ddim = make_indexdim(Dis,Dblockind);
            auto Dminm = std::numeric_limits<size_t>::max();
            for(decltype(Ddim.size()) j = 0; j < Ddim.size(); ++j)
                {
                Dminm = std::min(Dminm,Ddim[j]);
                }

            auto Dref = VecRefc(dblock,Dminm);

            contractDiagPartial(Dref,Dind,
                                Tref,Tind,
                                Cref,Cind);
            };

        loopContractedBlocks(D,Dis,
                             T,Tis,
                             C,Cis,
                             do_contract);
        }
    else
        {
        Error("Fully contracted IQTDiag*IQTReal not yet implemented");
        //auto nd = m.makeNewData<IQTDiag>(Cis,Cdiv);
        //auto& C = *nd;

        //loopContractedBlocks(D,Dis,
        //                     T,Tis,
        //                     C,Cis,
        //                     do_contract);
        }
    }

void
doTask(Contract<IQIndex>& Con,
       IQTDiag const& A,
       IQTReal const& B,
       ManageStore& m)
    {
    Label Aind,
          Bind,
          Cind;
    bool sortInds = false;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),Aind,Bind);
    contractIS(Con.Lis,Aind,Con.Ris,Bind,Con.Nis,Cind,sortInds);
    blockDiagDense(A,Con.Lis,Aind,
                   B,Con.Ris,Bind,
                   Con.Nis,Cind,m);
    }

void
doTask(Contract<IQIndex>& Con,
       IQTReal const& A,
       IQTDiag const& B,
       ManageStore& m)
    {
    Label Aind,
          Bind,
          Cind;
    bool sortInds = false;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),Aind,Bind);
    contractIS(Con.Lis,Aind,Con.Ris,Bind,Con.Nis,Cind,sortInds);
    blockDiagDense(B,Con.Ris,Bind,
                   A,Con.Lis,Aind,
                   Con.Nis,Cind,m);
    }

void
doTask(Conj, IQTDiag const& d) { }


Real
doTask(NormNoScale, 
       IQTDiag const& d) 
    { 
    Real nrm = 0;
    for(auto& elt : d.store)
        {
        nrm += elt*elt;
        }
    return std::sqrt(nrm);
    }


void
doTask(PrintIT<IQIndex> & P, 
       IQTDiag const& d)
    {
    P.s << "IQTDiag {" << d.offsets.size() << " blocks; Data size = " << d.store.size() << "}\n\n";
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
    for(auto& io : d.offsets)
        {
        bool indices_printed = false;
        //Determine block indices (where in the IQIndex space
        //this non-zero block is located)
        inverseBlockInd(io.block,P.is,block);
        //Wire up GCounter with appropriate dims
        auto blockm = blockIndex(0).m();
        for(auto i : count(1,rank)) blockm = std::min(blockm,blockIndex(i).m());

        auto os = io.offset;
        for(auto n : count(blockm))
            {
            auto val = scalefac*d.store[os++];
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
                for(auto ii = 1; ii <= rank; ++ii)
                    {
                    P.s << (1+n);
                    if(ii < rank) P.s << ",";
                    }
                P.s << ") ";
                P.printVal(val);
                }
            }
        }
    }

void
doTask(Write & W, IQTDiag const& d)
    {
    W.writeType(StorageType::IQTDiag,d); 
    }

} //namespace itensor

