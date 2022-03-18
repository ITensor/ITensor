//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//#include "itensor/util/iterate.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/algs.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"
#include "itensor/itdata/dense.h"
#include "itensor/itdata/qdense.h"
#include "itensor/itdata/qutil.h"
#include "itensor/util/print_macro.h"

using std::vector;
using std::string;
using std::move;

namespace itensor {

const char*
typeNameOf(QDenseReal const& d) { return "QDenseReal"; }
const char*
typeNameOf(QDenseCplx const& d) { return "QDenseCplx"; }

BlOf 
make_blof(Block const& b, long o)
    {
    BlOf B;
    B.block = b;
    B.offset = o;
    return B;
    }

bool
operator==(Block const& b1, Block const& b2)
    {
    for(auto i : range(b1.size()))
        {
        if(b1[i] != b2[i]) return false;
        }
    return true;
    }

bool
operator!=(Block const& b1, Block const& b2) { return !(b1==b2); }

bool
operator<(Block const& l1, Block const& l2)
    {
    return std::lexicographical_compare(l1.rbegin(),l1.rend(),
                                        l2.rbegin(),l2.rend());
    }

bool
operator>(Block const& l1, Block const& l2) { return !(l1 < l2) && (l1 != l2); }

//function object for calling binaryFind
//on offset vectors below
struct compBlock
    {
    bool
    operator()(const BlOf& bo1,
               const BlOf& bo2) const
        { return bo1.block < bo2.block; }
    bool
    operator()(const BlOf& bo, Block const& blk) const        
        { return bo.block < blk; }
    bool
    operator()(Block const& blk, const BlOf& bo) const 
        { return blk < bo.block; }
    };

QN
calcDiv(IndexSet const& is, 
        Block const& block_ind)
    {
    QN div;
    for(auto i : range(order(is))) { div += is[i].dir()*is[i].qn(1+block_ind[i]); }
    return div;
    }

template<typename T>
QN
doTask(CalcDiv const& C,
       QDense<T> const& D)
    {
    if(order(C.is)==0 || D.offsets.empty()) return QN{};
    auto b = D.offsets.front().block;
    auto block_ind = Block(order(C.is));
    block_ind = b;
    return calcDiv(C.is,block_ind);
    }
template QN doTask(CalcDiv const&,QDense<Real> const&);
template QN doTask(CalcDiv const&,QDense<Cplx> const&);

template<typename T>
QDense<T>::
QDense(IndexSet const& is, 
       QN       const& div)
    {
    auto totalsize = updateOffsets(is,div);
    store.assign(totalsize,0.);
    }
template QDense<Real>::QDense(IndexSet const&, QN const&);
template QDense<Cplx>::QDense(IndexSet const&, QN const&);

// Constructor taking a list of block labels
// instead of QN divergence
template<typename T>
QDense<T>::
QDense(IndexSet const& is,
       Blocks   const& blocks)
    {
    auto totalsize = updateOffsets(is,blocks);
    store.assign(totalsize,0.);
    }
template QDense<Real>::QDense(IndexSet const&, Blocks const&);
template QDense<Cplx>::QDense(IndexSet const&, Blocks const&);

std::tuple<BlockOffsets,long>
getBlockOffsets(IndexSet const& is,
                QN       const& div)
    {
    auto bofs = BlockOffsets();
    if(order(is)==0)
        {
        bofs.push_back(make_blof(Block(0),0));
        return std::make_tuple(bofs,1);
        }
    //Set up a Range to iterate over all blocks
    auto RB = RangeBuilder(order(is));
    for(auto j : range(order(is)))
        RB.nextIndex(is[j].nblock());

    long totalsize = 0;
    for(auto I : RB.build())
        {
        auto blockqn = QN{};
        for(auto j : range(order(is)))
            {
            auto& J = is[j];
            blockqn += J.qn(1+I[j])*J.dir();
            }
        if(blockqn == div)
            {
            auto block = Block(order(is));
            long totdim = 1;   //accumulate dim of Indices
            for(auto j : range(order(is)))
                {
                auto& J = is[j];
                auto i_j = I[j];
                block[j] = i_j;
                totdim *= J.blocksize0(i_j);
                }
            bofs.push_back(make_blof(block,totalsize));
            totalsize += totdim;
            }
        }
    return std::make_tuple(bofs,totalsize);
    }

template<typename T>
long QDense<T>::
updateOffsets(IndexSet const& is,
              QN       const& div)
    {
    auto [bofs,size] = getBlockOffsets(is,div);
    offsets = bofs;
    return size;
    }

template<typename T>
long QDense<T>::
updateOffsets(IndexSet const& is,
              Blocks   const& blocks)
    {
    offsets.clear();

    if(order(is)==0)
        {
        offsets.push_back(make_blof(Block(0),0));
        return 1;
        }

    long totalsize = 0;
    for(auto const& block : blocks)
        {
        long totdim = 1;   //accumulate dim of Indices
        for(auto j : range(order(is)))
            {
            auto& J = is[j];
            auto i_j = block[j];
            totdim *= J.blocksize0(i_j);
            }
        offsets.push_back(make_blof(block,totalsize));
        totalsize += totdim;
        }
    return totalsize;
    }

long
offsetOf(BlockOffsets const& offsets,
         Block        const& blockind)
    {
    auto blk = detail::binaryFind(offsets,blockind,compBlock());
    if(blk) return blk->offset;
    return -1;
    }

int
offsetOfLoc(BlockOffsets const& offsets,
            Block        const& blockind)
    {
    auto it = std::lower_bound(offsets.begin(),offsets.end(),blockind,compBlock());
    int loc = std::distance(offsets.begin(),it);
    return loc;
    }

Cplx
doTask(GetElt& G, QDenseReal const& d)
    {
    auto* pelt = d.getElt(G.is,G.inds);
    if(pelt) return Cplx(*pelt,0.);
    return Cplx(0.,0.);
    }
Cplx
doTask(GetElt& G, QDenseCplx const& d)
    {
    auto* pelt = d.getElt(G.is,G.inds);
    if(pelt) return *pelt;
    return Cplx(0.,0.);
    }

template<typename E, typename T>
void
setEltImpl(SetElt<E> & S, QDense<T> & d)
    {
    auto eltblockoffset = d.getEltBlockOffset(S.is,S.inds);
    auto* pelt = std::get<0>(eltblockoffset);
    auto block = std::get<1>(eltblockoffset);
    auto eltoffset = std::get<2>(eltblockoffset);
    if(pelt)
      {
      // The block already exists
      *pelt = S.elt;
      }
    else
      {
      // The block doesn't exist, so add it and then
      // set the element
      auto boffset = d.insertBlock(S.is,block);
      pelt = d.store.data()+boffset+eltoffset;
      *pelt = S.elt;
      }
    }

template<typename T>
void
doTask(SetElt<Real>& S, QDense<T>& d)
    {
    setEltImpl<Real,T>(S,d);
    }
template void doTask(SetElt<Real>&, QDense<Real>&);
template void doTask(SetElt<Real>&, QDense<Cplx>&);

void
doTask(SetElt<Cplx>& S, QDenseCplx & d)
    {
    setEltImpl<Cplx,Cplx>(S,d);
    }

void
doTask(SetElt<Cplx>& S, QDenseReal const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDenseCplx>(d.offsets,d.begin(),d.end());
    setEltImpl<Cplx,Cplx>(S,*nd);
    }

template<typename T>
Cplx
doTask(SumEls, QDense<T> const& d)
    {
    Cplx s = 0.;
    for(auto& el : d.store) s += el;
    return s;
    }
template Cplx doTask(SumEls, QDense<Real> const&);
template Cplx doTask(SumEls, QDense<Cplx> const&);

template<typename T>
void
doTask(Mult<Real> const& M, QDense<T>& D)
    {
    auto d = realData(D);
    dscal_wrapper(d.size(),M.x,d.data());
    }
template void doTask(Mult<Real> const&, QDenseReal&);
template void doTask(Mult<Real> const&, QDenseCplx&);


void
doTask(Mult<Cplx> const& M, QDense<Cplx> & d)
    {
    for(auto& el : d) el *= M.x;
    }

void
doTask(Mult<Cplx> const& M, QDense<Real> const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDenseCplx>(d.offsets,d.begin(),d.end());
    doTask(M,*nd);
    }

void
doTask(MakeCplx const&, QDense<Cplx> & d)
    {
    //nothing to do
    }
void
doTask(MakeCplx const&, QDense<Real> const& d, ManageStore & m)
    {
    m.makeNewData<QDenseCplx>(d.offsets,d.begin(),d.end());
    }

template<typename T>
int
doTask(NNZBlocks, QDense<T> const& D)
    {
    return D.offsets.size();
    }
template int doTask(NNZBlocks, QDense<Real> const&);
template int doTask(NNZBlocks, QDense<Cplx> const&);

template<typename T>
long
doTask(NNZ, QDense<T> const& D)
    {
    return D.store.size();
    }
template long doTask(NNZ, QDense<Real> const&);
template long doTask(NNZ, QDense<Cplx> const&);

template<typename T>
void
doTask(Fill<T> const& F, QDense<T> & d)
    {
    stdx::fill(d,F.x);
    }
template void doTask(Fill<Real> const&, QDense<Real> &);
template void doTask(Fill<Cplx> const&, QDense<Cplx> &);

template<typename FT, typename DT,class>
void
doTask(Fill<FT> const& F, QDense<DT> const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDense<FT>>(d.offsets,d.size());
    doTask(F,*nd);
    }
template void doTask(Fill<Real> const& F, QDense<Cplx> const&, ManageStore &);
template void doTask(Fill<Cplx> const& F, QDense<Real> const&, ManageStore &);


void
doTask(Conj, QDenseCplx & d)
    {
    for(auto& el : d) applyConj(el);
    }

void
doTask(TakeReal, QDenseCplx const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDenseReal>(d.offsets,d.size());
    for(auto i : range(d))
        {
        nd->store[i] = d.store[i].real();
        }
    }

void
doTask(TakeImag, QDenseReal & d)
    { 
    //Set all elements to zero
    doTask(Fill<Real>{0.},d);
    }

void
doTask(TakeImag, QDenseCplx const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDenseReal>(d.offsets,d.size());
    for(auto i : range(d))
        {
        nd->store[i] = d.store[i].imag();
        }
    }

template<typename T>
Real
doTask(NormNoScale, QDense<T> const& D)
    { 
    auto d = realData(D);
    return dnrm2_wrapper(d.size(),d.data());
    }
template Real doTask(NormNoScale, QDense<Real> const& D);
template Real doTask(NormNoScale, QDense<Cplx> const& D);

template<typename T>
void
doTask(PrintIT& P, QDense<T> const& d)
    {
    auto name = format("QDense %s",typeName<T>());
    if(not P.print_data)
        {
        P.printInfo(d,name,doTask(NormNoScale{},d));
        return;
        }

    P.s << format("QDense %s {%d blocks; data size %d}\n",
                  typeName<T>(),d.offsets.size(),d.size());
    //Real scalefac = 1.0;
    //if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    //else P.s << "(omitting too large scale factor)\n";

    auto ord = order(P.is);
    if(ord == 0) 
        {
        P.s << "  ";
        P.s << formatVal(d.store.front()) << "\n";
        return;
        }
        
    Range brange;
    auto C = detail::GCounter(ord);
    for(auto const& io : d.offsets)
        {
        bool block_info_printed = false;

        Block boff(ord,0);
        for(auto i : range(ord))
            {
            for(auto j : range(io.block[i]))
                boff[i] += P.is[i].blocksize0(j);
            }

        //Wire up GCounter with appropriate dims
        C.reset();
        for(auto i : range(ord))
            {
            C.setRange(i,0,P.is[i].blocksize0(io.block[i])-1);
            }
        for(auto os = io.offset; C.notDone(); ++C, ++os)
            {
            auto val = d.store[os];
            if(std::norm(val) >= Global::printScale())
                {
                if(not block_info_printed)
                    {
                    block_info_printed = true;
                    //Print Indices of this block
                    P.s << "Block:";
                    for(auto const& bi : io.block)
                        {
                        P.s << " " << (1+bi);
                        //if(i > 0) P.s << " ";
                        //P.s << blockIndex(i) << "<" << P.is[i].dir() << ">";
                        }
                    P.s << "\n";
                    }

                P.s << "(";
                for(auto ii : range(ord))
                    {
                    P.s << (1+boff[ii]+C[ii]);
                    if(1+ii != ord) P.s << ",";
                    }
                P.s << ") ";

                //P.s << "[";
                //for(auto ii : range(ord))
                //    {
                //    P.s << (1+C[ii]);
                //    if(1+ii != ord) P.s << ",";
                //    }
                //P.s << "] ";

                P.s << formatVal(val) << "\n";
                }
            }
        }
    }
template void doTask(PrintIT& P, QDense<Real> const& d);
template void doTask(PrintIT& P, QDense<Cplx> const& d);

struct Adder
    {
    const Real f = 1.;
    Adder(Real f_) : f(f_) { }
    template<typename T1, typename T2>
    void operator()(T2 v2, T1& v1) { v1 += f*v2; }
    void operator()(Cplx v2, Real& v1) { }
    };

template<typename T1, typename T2>
void
add(PlusEQ const& P,
    QDense<T1>            & A,
    QDense<T2>       const& B)
    {
    auto r = order(P.is1());

    if(r==0)
        {
        auto dA = realData(A);
        auto dB = realData(B);
        daxpy_wrapper(dA.size(),P.alpha(),dB.data(),1,dA.data(),1);
        return;
        }

    auto Bblock = Block(r,0);
    Range Arange,
          Brange;

    for(auto const& aio : A.offsets)
        {
        for(auto i : range(r))
            Bblock[i] = aio.block[P.perm().dest(i)];

        auto bblock = getBlock(B,P.is2(),Bblock);
        if(!bblock) continue;

        Arange.init(make_indexdim(P.is1(),aio.block));
        Brange.init(make_indexdim(P.is2(),Bblock));
        auto aref = makeTenRef(A.data(),aio.offset,A.size(),&Arange);
        auto bref = makeRef(bblock,&Brange);
        transform(permute(bref,P.perm()),aref,Adder{P.alpha()});
        }
    }

template<typename TA, typename TB>
void
doTask(PlusEQ const& P,
       QDense<TA>      const& A,
       QDense<TB>      const& B,
       ManageStore          & m)
    {
    if(B.store.size() == 0) return;

    //
    // If B has blocks that A doesn't have,
    // then we need to widen the storage of A
    //

    auto r = order(P.is1());

    if(r == 0)
        {
        if(isReal(A) && isCplx(B))
            {
            auto *nA = m.makeNewData<QDenseCplx>(A.offsets,A.begin(),A.end());
            add(P,*nA,B);
            }
        else
            {
            auto *mA = m.modifyData(A);
            add(P,*mA,B);
            }
        return;
        }

    // Store the blocks of the output
    // TODO: can this be optimized more? Currently, the
    // strategy is to permute and sort the blocks of B,
    // then mergy A and B using that they are both sorted
    auto Cblocks = Blocks();
    // Reserve the maximum space we need in order to avoid
    // reallocations when using push_back()
    Cblocks.reserve(A.offsets.size()+B.offsets.size());

    // First we need to permute the blocks of B
    // and sort them
    auto Bblockps = Blocks(B.offsets.size(),Block(r));
    auto invperm = inverse(P.perm());
    for(auto ib : range(B.offsets.size()))
        {
        auto const& Bblock = B.offsets[ib].block;
        auto& Bblockp = Bblockps[ib];
        for(auto i : range(r))
            Bblockp[i] = Bblock[invperm.dest(i)];
        }
    std::sort(Bblockps.begin(),Bblockps.end());

    size_t ia = 0,
           ib = 0;
    while(ia < A.offsets.size() && ib < B.offsets.size())
        {
        auto const& Ablock = A.offsets[ia].block;
        auto const& Bblockp = Bblockps[ib];
        if(Bblockp < Ablock)
            {
            Cblocks.push_back(Bblockp);
            ib++;
            }
        else if(Ablock < Bblockp)
            {
            Cblocks.push_back(Ablock);
            ia++;
            }
        else // Ablock == Bblockp
            {
            Cblocks.push_back(Ablock);
            ia++;
            ib++;
            }
        }
    while(ia < A.offsets.size())
        {
        auto const& Ablock = A.offsets[ia].block;
        Cblocks.push_back(Ablock);
        ia++;
        }
    while(ib < B.offsets.size())
        {
        auto const& Bblockp = Bblockps[ib];
        Cblocks.push_back(Bblockp);
        ib++;
        }

    // TODO: make a special case for B.offsets.size() == Cblocks.size()?
    //       This could avoid having to allocate new memory in certain
    //       situations
    if(A.offsets.size() < Cblocks.size())
        {
        // This means there are blocks in B that are not in A
        // Need to expand the data
        auto *nA = m.makeNewData<QDense<common_type<TA,TB>>>(P.is1(),Cblocks);
        // Do a trivial permutation
        auto trivial_perm = PlusEQ::permutation(r);
        auto PA = PlusEQ(trivial_perm,P.is1(),P.is1(),1.0);
        add(PA,*nA,A);
        add(P,*nA,B);
        }
    else if(isReal(A) && isCplx(B))
        {
        auto *nA = m.makeNewData<QDenseCplx>(A.offsets,A.begin(),A.end());
        add(P,*nA,B);
        }
    else
        {
        auto *mA = m.modifyData(A);
        add(P,*mA,B);
        }
    }
template void doTask(PlusEQ const&, QDense<Real> const&, QDense<Real> const&, ManageStore&);
template void doTask(PlusEQ const&, QDense<Real> const&, QDense<Cplx> const&, ManageStore&);
template void doTask(PlusEQ const&, QDense<Cplx> const&, QDense<Real> const&, ManageStore&);
template void doTask(PlusEQ const&, QDense<Cplx> const&, QDense<Cplx> const&, ManageStore&);


template<typename VA, typename VB>
void
doTask(Contract& Con,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m)
    {
    using VC = common_type<VA,VB>;
    Labels Lind,
           Rind;

    computeLabels(Con.Lis,order(Con.Lis),Con.Ris,order(Con.Ris),Lind,Rind);

    //compute new index set (Con.Nis):
    Labels Cind;
    const bool sortResult = false;
    contractIS(Con.Lis,Lind,Con.Ris,Rind,Con.Nis,Cind,sortResult);

    //Allocate storage for C
    auto [Coffsets,Csize,blockContractions] = getContractedOffsets(A,Con.Lis,B,Con.Ris,Con.Nis);

    // Create QDense storage with uninitialized memory, faster than
    // setting to zeros
    auto nd = m.makeNewData<QDense<VC>>(undef,Coffsets,Csize);
    auto& C = *nd;

    //Determines if the contraction in the list overwrites or
    //adds to the data. Initially, overwrite the data since the
    //data starts uninitialized
    auto betas = std::vector<Real>(C.offsets.size(),0.);

    //Function to execute for each pair of
    //contracted blocks of A and B
    auto do_contract = 
        [&Con,&Lind,&Rind,&Cind,&betas]
        (DataRange<const VA> ablock, Block const& Ablockind,
         DataRange<const VB> bblock, Block const& Bblockind,
         DataRange<VC>       cblock, Block const& Cblockind,
         int Cblockloc)
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
        auto aref = makeRef(ablock,&Arange);
        auto bref = makeRef(bblock,&Brange);
        auto cref = makeRef(cblock,&Crange);

        // cref += aref*bref or cref = aref*bref
        contract(aref,Lind,bref,Rind,cref,Cind,1.,betas[Cblockloc]);

        // If the block had not been called, betas[Cblockloc] == 0
        // Set it to 1 after it has been called
        betas[Cblockloc] = 1.;
        };

    loopContractedBlocks(A,Con.Lis,
                         B,Con.Ris,
                         C,Con.Nis,
                         blockContractions,
                         do_contract);

#ifdef USESCALE
    Con.scalefac = computeScalefac(C);
#endif
    }
template void doTask(Contract& Con,QDense<Real> const&,QDense<Real> const&,ManageStore&);
template void doTask(Contract& Con,QDense<Cplx> const&,QDense<Real> const&,ManageStore&);
template void doTask(Contract& Con,QDense<Real> const&,QDense<Cplx> const&,ManageStore&);
template void doTask(Contract& Con,QDense<Cplx> const&,QDense<Cplx> const&,ManageStore&);

template<typename VA, typename VB>
void
doTask(NCProd& P,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m)
    {
    using VC = common_type<VA,VB>;
    auto& Ais = P.Lis;
    auto& Bis = P.Ris;
    auto& Cis = P.Nis;
    auto rA = order(Ais);
    auto rB = order(Bis);
    Labels Aind,
           Bind,
           Cind;
    computeLabels(Ais,rA,Bis,rB,Aind,Bind);
    ncprod(Ais,Aind,Bis,Bind,Cis,Cind);

    Labels BtoA(rA,-1);
    for(auto ia : range(rA))
    for(auto ib : range(rB))
        if(Bis[ib] == Ais[ia])
            {
            BtoA[ib] = ia;
            break;
            }

    auto Cdiv = QN{};
        {
        Cdiv = doTask(CalcDiv{Ais},A);
        auto Ablock_ind = Block(rA);
        Ablock_ind = A.offsets.front().block;
        auto Bblock_ind = Block(rB);
        for(auto& bo : B.offsets)
            {
            Bblock_ind = bo.block;
            bool matchesA = true;
            for(auto n : range(rB))
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
        for(auto n : range(rB))
            if(Bind[n] > 0) //unique
                {
                Cdiv += Bis[n].dir()*Bis[n].qn(1+Bblock_ind[n]);
                }
        }

    //Allocate storage for C
    auto& C = *m.makeNewData<QDense<VC>>(Cis,Cdiv);

    auto do_ncprod = 
        [&P,&Aind,&Bind,&Cind]
        (DataRange<const VA> ablock, Block const& Ablockind,
         DataRange<const VB> bblock, Block const& Bblockind,
         DataRange<VC>       cblock, Block const& Cblockind)
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
        auto aref = makeRef(ablock,&Arange);
        auto bref = makeRef(bblock,&Brange);
        auto cref = makeRef(cblock,&Crange);

        //Compute cref += aref*bref
        ncprod(aref,Aind,bref,Bind,cref,Cind);
        };

    loopContractedBlocks(A,Ais,
                         B,Bis,
                         C,Cis,
                         do_ncprod);

#ifdef USESCALE
    P.scalefac = computeScalefac(C);
#endif
    }
template void doTask(NCProd&,QDense<Real> const&,QDense<Real> const&,ManageStore&);
template void doTask(NCProd&,QDense<Cplx> const&,QDense<Real> const&,ManageStore&);
template void doTask(NCProd&,QDense<Real> const&,QDense<Cplx> const&,ManageStore&);
template void doTask(NCProd&,QDense<Cplx> const&,QDense<Cplx> const&,ManageStore&);

template<typename T>
void
permuteQDense(Permutation  const& P,
              QDense<T>    const& dA,
              IndexSet   const& Ais,
              QDense<T>         & dB,
              IndexSet        & Bis)
    {
    // Recalculate new indexset by permuting
    // original indexset (otherwise it segfaults)
    auto r = order(Ais);
    auto bind = IndexSetBuilder(r);
    for(auto i : range(r))
        bind.setIndex(P.dest(i),Ais[i]);
    Bis = bind.build();
    dB = QDense<T>(Bis,doTask(CalcDiv{Ais},dA));
    // Perform permutation
    auto Bblock = Block(r,-1);
    Range Arange,
          Brange;
    for(auto const& aio : dA.offsets)
        {
        //Compute bi, new block index of blk
        for(auto j : range(aio.block))
            Bblock.at(P.dest(j)) = aio.block[j];
        Arange.init(make_indexdim(Ais,aio.block));
        Brange.init(make_indexdim(Bis,Bblock));

        auto bblock = getBlock(dB,Bis,Bblock);
        auto bref = makeRef(bblock,&Brange);
        auto aref = makeTenRef(dA.data(),aio.offset,dA.size(),&Arange);

        bref += permute(aref,P);
        }
    }

template<typename T>
void
doTask(Order const& O,
       QDense<T> & dB)
    {
    auto const dA = dB;
    auto Bis = O.is2();
    permuteQDense(O.perm(),dA,O.is1(),dB,Bis);
    }
template void doTask(Order const&,QDense<Real> &);
template void doTask(Order const&,QDense<Cplx> &);

template<typename V>
TenRef<Range,V>
doTask(GetBlock<V> const& G,
       QDense<V> & d)
    {
    auto block = getBlock(d,G.is,G.block_ind);
    auto RB = RangeBuilder(order(G.is));
    for(auto j : range(order(G.is)))
        {
        RB.nextIndex(G.is[j].blocksize0(G.block_ind[j]));
        }
    return makeRef(block,RB.build());
    }
template TenRef<Range,Real> doTask(GetBlock<Real> const& G,QDense<Real> & d);
template TenRef<Range,Cplx> doTask(GetBlock<Cplx> const& G,QDense<Cplx> & d);

template<typename T>
bool
doTask(IsDense,
       QDense<T> const& d)
    {
    return true;
    }
template bool doTask(IsDense,QDense<Real> const& d);
template bool doTask(IsDense,QDense<Cplx> const& d);

template<typename V>
void
doTask(RemoveQNs & R, 
       QDense<V> const& d,
       ManageStore & m)
    {
    auto r = order(R.is);
    auto *nd = m.makeNewData<Dense<V>>(dim(R.is),0);
    auto *pd = d.data();
    auto *pn = nd->data();
    detail::GCounter C(r);
    for(auto const& io : d.offsets)
        {
        for(auto j : range(r))
            {
            long start = 0;
            for(auto b : range(io.block[j]))
                {
                start += R.is[j].blocksize0(b);
                }
            C.setRange(j,start,start+R.is[j].blocksize0(io.block[j])-1);
            }
        //TODO: need to make a Range/TensorRef iterator
        //to rewrite the following code more efficiently
        for(; C.notDone(); ++C)
            {
            pn[offset(R.is,C.i)] = pd[io.offset+C.ind];
            }
        }
    }
template void doTask(RemoveQNs &, QDense<Real> const&, ManageStore &);
template void doTask(RemoveQNs &, QDense<Cplx> const&, ManageStore &);

std::ostream&
operator<<(std::ostream & s, BlOf const& blof)
    {
    s << "Block: " << blof.block << ", Offset: " << blof.offset << "\n";
    return s;
    }

std::ostream&
operator<<(std::ostream & s, BlockOffsets const& offsets)
    {
    for(auto const& blof : offsets) s << blof;
    return s;
    }

std::ostream&
operator<<(std::ostream & s, Blocks const& blocks)
    {
    for(auto const& block : blocks) s << block << "\n";
    return s;
    }

template<typename T>
std::ostream&
operator<<(std::ostream & s, QDense<T> const& t)
    {
    s << "QDense blocks and offsets:\n";
    s << t.offsets << "\n";
    s << "\nQDense storage:\n";
    for(auto i : range(t.store.size()))
      s << i << " " << t.store[i] << "\n";
    return s;
    }
template std::ostream& operator<<(std::ostream & s, QDense<Real> const& t);
template std::ostream& operator<<(std::ostream & s, QDense<Cplx> const& t);

#ifdef ITENSOR_USE_HDF5

vector<long>
offsets_to_array(BlockOffsets const& boff, int N)
    {
    auto nblocks = boff.size();
    auto asize = (N+1)*nblocks;
    auto n = 0;
    auto a = vector<long>(asize);
    for(auto& bo : boff)
        {
        for(auto j : range(N))
            {
            a[n] = bo.block[j]+1;
            n += 1;
            }
        a[n] = bo.offset;
        n += 1;
        }
    return a;
    }

BlockOffsets
array_to_offsets(vector<long> const& a, long N)
    {
    auto asize = a.size();
    auto nblocks = ldiv(asize, N+1).quot;
    auto boff = BlockOffsets(nblocks);
    long j = 0;
    for(auto n : range(nblocks))
        {
        auto block = Block(N);
        for(auto m : range(N)) block[m] = a[j+m]-1;
        long offset = a[j+N];
        boff[n] = BlOf{block,offset};
        j += (N + 1);
        }
    return boff;
    }

const char*
juliaTypeNameOf(QDenseReal const& d) { return "BlockSparse{Float64}"; }
const char*
juliaTypeNameOf(QDenseCplx const& d) { return "BlockSparse{ComplexF64}"; }

template<typename V>
void
h5_write(h5::group parent, std::string const& name, QDense<V> const& D)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type",juliaTypeNameOf(D),true);
    h5_write_attribute(g,"version",long(1));
    long N = 0;
    if(!D.offsets.empty()) N = D.offsets.front().block.size();
    h5_write(g,"ndims",N);
    auto off_array = offsets_to_array(D.offsets,N);
    h5_write(g,"offsets",off_array);
    auto data = std::vector<V>(D.store.begin(),D.store.end());
    h5_write(g,"data",data);
    }
template void h5_write(h5::group, std::string const&, QDense<Real> const& D);
template void h5_write(h5::group, std::string const&, QDense<Cplx> const& D);

template<typename V>
void
h5_read(h5::group parent, std::string const& name, QDense<V> & D)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<string>(g,"type");
    if(type != juliaTypeNameOf(D)) 
        {
        Error(format("Group does not contain %s or %s data in HDF5 file",typeNameOf(D),juliaTypeNameOf(D)));
        }
    auto N = h5_read<long>(g,"ndims");
    auto off_array = offsets_to_array(D.offsets,N);
    auto offsets = h5_read<vector<long>>(g,"offsets");
    auto boff = array_to_offsets(offsets,N);
    auto data = h5_read<vector<V>>(g,"data");
    D = QDense(boff,data);
    }
template void h5_read(h5::group, std::string const&, QDense<Real> & D);
template void h5_read(h5::group, std::string const&, QDense<Cplx> & D);


#endif //ITENSOR_USE_HDF5

} //namespace itensor

