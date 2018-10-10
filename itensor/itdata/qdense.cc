//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
//#include "itensor/util/range.h"
#include "itensor/detail/gcounter.h"
#include "itensor/detail/algs.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/tensor/contract.h"
#include "itensor/itdata/qdense.h"
#include "itensor/itdata/qutil.h"
#include "itensor/util/print_macro.h"

using std::vector;
using std::move;

namespace itensor {

BlOf inline
make_blof(long b, long o)
    {
    BlOf B;
    B.block = b;
    B.offset = o;
    return B;
    }

const char*
typeNameOf(QDenseReal const& d) { return "QDenseReal"; }
const char*
typeNameOf(QDenseCplx const& d) { return "QDenseCplx"; }

//function object for calling binaryFind
//on offset vectors below
struct compBlock
    {
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
        Labels const& block_ind)
    {
    QN div;
    for(auto i : range(is.r())) { div += is[i].dir()*is[i].qn(1+block_ind[i]); }
    return div;
    }

template<typename T>
QN
doTask(CalcDiv const& C,
       QDense<T> const& D)
    {
    if(C.is.r()==0 || D.offsets.empty()) return QN{};
    auto b = D.offsets.front().block;
    Labels block_ind(C.is.r());
    computeBlockInd(b,C.is,block_ind);
    return calcDiv(C.is,block_ind);
    }
template QN doTask(CalcDiv const&,QDense<Real> const&);
template QN doTask(CalcDiv const&,QDense<Cplx> const&);

template<typename T>
QDense<T>::
QDense(IQIndexSet const& is, 
       QN         const& div)
    {
    auto totalsize = updateOffsets(is,div);
    store.assign(totalsize,0.);
    }
template QDense<Real>::QDense(IQIndexSet const&, QN const&);
template QDense<Cplx>::QDense(IQIndexSet const&, QN const&);

template<typename T>
long QDense<T>::
updateOffsets(IQIndexSet const& is,
              QN         const& div)
    {
    offsets.clear();

    if(is.r()==0)
        {
        offsets.push_back(make_blof(0,0));
        return 1;
        }

    //Set up a Range to iterate over all blocks
    auto RB = RangeBuilder(is.r());
    for(auto j : range(is.r()))
        RB.nextIndex(is[j].nindex());

    long totalsize = 0;
    for(auto I : RB.build())
        {
        auto blockqn = QN{};
        for(auto j : range(is.r()))
            {
            auto& J = is[j];
            blockqn += J.qn(1+I[j])*J.dir();
            }
        if(blockqn == div)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 totm = 1;   //accumulate area of Indices
            for(auto j : range(is.r()))
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
offsetOf(std::vector<BlOf> const& offsets,
         long blockind)
    {
    auto blk = detail::binaryFind(offsets,blockind,compBlock());
    if(blk) return blk->offset;
    return -1;
    }

Cplx
doTask(GetElt<IQIndex>& G, QDenseReal const& d)
    {
    auto* pelt = d.getElt(G.is,G.inds);
    if(pelt) return Cplx(*pelt,0.);
    return Cplx(0.,0.);
    }
Cplx
doTask(GetElt<IQIndex>& G, QDenseCplx const& d)
    {
    auto* pelt = d.getElt(G.is,G.inds);
    if(pelt) return *pelt;
    return Cplx(0.,0.);
    }

template<typename E, typename T>
void
setEltImpl(SetElt<E,IQIndex> & S, QDense<T> & d)
    {
    auto* pelt = d.getElt(S.is,S.inds);
    if(pelt) *pelt = S.elt;
    else     Error("Setting IQTensor element non-zero would violate its symmetry.");
    }

template<typename T>
void
doTask(SetElt<Real,IQIndex>& S, QDense<T>& d)
    {
    setEltImpl<Real,T>(S,d);
    }
template void doTask(SetElt<Real,IQIndex>&, QDense<Real>&);
template void doTask(SetElt<Real,IQIndex>&, QDense<Cplx>&);

void
doTask(SetElt<Cplx,IQIndex>& S, QDenseCplx & d)
    {
    setEltImpl<Cplx,Cplx>(S,d);
    }

void
doTask(SetElt<Cplx,IQIndex>& S, QDenseReal const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDenseCplx>(d.offsets,d.begin(),d.end());
    setEltImpl<Cplx,Cplx>(S,*nd);
    }

template<typename T>
Cplx
doTask(SumEls<IQIndex>, QDense<T> const& d)
    {
    Cplx s = 0.;
    for(auto& el : d.store) s += el;
    return s;
    }
template Cplx doTask(SumEls<IQIndex>, QDense<Real> const&);
template Cplx doTask(SumEls<IQIndex>, QDense<Cplx> const&);

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
doTask(PrintIT<IQIndex>& P, QDense<T> const& d)
    {
    P.s << format("QDense %s {%d blocks; data size %d}\n",
                  typeName<T>(),d.offsets.size(),d.size());
    Real scalefac = 1.0;
    if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    else P.s << "(omitting too large scale factor)\n";

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        P.s << formatVal(scalefac*d.store.front()) << "\n";
        return;
        }
        
    Labels block(rank,0);
    auto blockIndex = [&block,&P](long i)->Index { return (P.is[i])[block[i]]; };

    Range brange;
    auto C = detail::GCounter(rank);
    for(const auto& io : d.offsets)
        {
        bool indices_printed = false;
        //Determine block indices (where in the IQIndex space
        //this non-zero block is located)
        computeBlockInd(io.block,P.is,block);

        Labels boff(rank,0);
        for(auto i : range(rank))
            {
            for(auto j : range(block[i]))
                boff[i] += P.is[i][j].m();
            }

        //Wire up GCounter with appropriate dims
        C.reset();
        for(decltype(rank) i = 0; i < rank; ++i)
            C.setRange(i,0,blockIndex(i).m()-1);
        for(auto os = io.offset; C.notDone(); ++C, ++os)
            {
            auto val = scalefac*d.store[os];
            if(std::norm(val) >= Global::printScale())
                {
                if(!indices_printed)
                    {
                    indices_printed = true;
                    //Print Indices of this block
                    for(auto i : range(rank))
                        {
                        if(i > 0) P.s << " ";
                        P.s << blockIndex(i) << "<" << P.is[i].dir() << ">";
                        }
                    P.s << "\n";
                    }
                P.s << "(";
                for(auto ii : range(rank))
                    {
                    P.s << (1+boff[ii]+C[ii]);
                    if(1+ii != rank) P.s << ",";
                    }
                P.s << ") ";

                //P.s << "[";
                //for(auto ii : range(rank))
                //    {
                //    P.s << (1+C[ii]);
                //    if(1+ii != rank) P.s << ",";
                //    }
                //P.s << "] ";

                P.s << formatVal(val) << "\n";
                }
            }
        }
    }
template void doTask(PrintIT<IQIndex>& P, QDense<Real> const& d);
template void doTask(PrintIT<IQIndex>& P, QDense<Cplx> const& d);

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
add(PlusEQ<IQIndex> const& P,
    QDense<T1>            & A,
    QDense<T2>       const& B)
    {
#ifdef DEBUG
    if(A.store.size() != B.store.size()) Error("Mismatched sizes in plusEq");
#endif
    if(isTrivial(P.perm()) && std::is_same<T1,T2>::value)
        {
        auto dA = realData(A);
        auto dB = realData(B);
        daxpy_wrapper(dA.size(),P.alpha(),dB.data(),1,dA.data(),1);
        }
    else
        {
        auto r = P.is1().r();
        Labels Ablock(r,0),
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
            auto bref = makeRef(bblock,&Brange);
            transform(permute(bref,P.perm()),aref,Adder{P.alpha()});
            }
        }
    }

template<typename TA, typename TB>
void
doTask(PlusEQ<IQIndex> const& P,
       QDense<TA>      const& A,
       QDense<TB>      const& B,
       ManageStore          & m)
    {
    if(B.store.size() == 0) return;

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
    }
template void doTask(PlusEQ<IQIndex> const&, QDense<Real> const&, QDense<Real> const&, ManageStore&);
template void doTask(PlusEQ<IQIndex> const&, QDense<Real> const&, QDense<Cplx> const&, ManageStore&);
template void doTask(PlusEQ<IQIndex> const&, QDense<Cplx> const&, QDense<Real> const&, ManageStore&);
template void doTask(PlusEQ<IQIndex> const&, QDense<Cplx> const&, QDense<Cplx> const&, ManageStore&);


template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m)
    {
    using VC = common_type<VA,VB>;
    Labels Lind,
          Rind;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),Lind,Rind);
    //compute new index set (Con.Nis):
    Labels Cind;
    const bool sortResult = false;
    contractIS(Con.Lis,Lind,Con.Ris,Rind,Con.Nis,Cind,sortResult);

    auto Cdiv = doTask(CalcDiv{Con.Lis},A)+doTask(CalcDiv{Con.Ris},B);

    //Allocate storage for C
    START_TIMER(33)
    auto nd = m.makeNewData<QDense<VC>>(Con.Nis,Cdiv);
    STOP_TIMER(33)
    auto& C = *nd;

    //Function to execute for each pair of
    //contracted blocks of A and B
    auto do_contract = 
        [&Con,&Lind,&Rind,&Cind]
        (DataRange<const VA> ablock, Labels const& Ablockind,
         DataRange<const VB> bblock, Labels const& Bblockind,
         DataRange<VC>       cblock, Labels const& Cblockind)
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

#ifdef USESCALE
    START_TIMER(21)
    Con.scalefac = computeScalefac(C);
    STOP_TIMER(21)
#endif
    }
template void doTask(Contract<IQIndex>& Con,QDense<Real> const&,QDense<Real> const&,ManageStore&);
template void doTask(Contract<IQIndex>& Con,QDense<Cplx> const&,QDense<Real> const&,ManageStore&);
template void doTask(Contract<IQIndex>& Con,QDense<Real> const&,QDense<Cplx> const&,ManageStore&);
template void doTask(Contract<IQIndex>& Con,QDense<Cplx> const&,QDense<Cplx> const&,ManageStore&);

template<typename VA, typename VB>
void
doTask(NCProd<IQIndex>& P,
       QDense<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m)
    {
    using VC = common_type<VA,VB>;
    auto& Ais = P.Lis;
    auto& Bis = P.Ris;
    auto& Cis = P.Nis;
    auto rA = rank(Ais);
    auto rB = rank(Bis);
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
        auto Ablock_ind = Labels(rA);
        computeBlockInd(A.offsets.front().block,Ais,Ablock_ind);
        auto Bblock_ind = Labels(rB);
        for(auto& bo : B.offsets)
            {
            computeBlockInd(bo.block,Bis,Bblock_ind);
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
        (DataRange<const VA> ablock, Labels const& Ablockind,
         DataRange<const VB> bblock, Labels const& Bblockind,
         DataRange<VC>       cblock, Labels const& Cblockind)
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
template void doTask(NCProd<IQIndex>&,QDense<Real> const&,QDense<Real> const&,ManageStore&);
template void doTask(NCProd<IQIndex>&,QDense<Cplx> const&,QDense<Real> const&,ManageStore&);
template void doTask(NCProd<IQIndex>&,QDense<Real> const&,QDense<Cplx> const&,ManageStore&);
template void doTask(NCProd<IQIndex>&,QDense<Cplx> const&,QDense<Cplx> const&,ManageStore&);

template<typename T>
void
permuteQDense(Permutation  const& P,
              QDense<T>    const& dA,
              IQIndexSet   const& Ais,
              QDense<T>         & dB,
              IQIndexSet        & Bis)
    {
    // Recalculate new indexset by permuting
    // original indexset (otherwise it segfaults)
    auto r = Ais.r();
    auto bind = IQIndexSetBuilder(r);
    for(auto i : range(r))
        {
        bind.setIndex(P.dest(i),Ais[i]);
        }
    Bis = bind.build();
    dB = QDense<T>(Bis,doTask(CalcDiv{Ais},dA));
    // Perform permutation
    Labels Ablock(r,-1),
          Bblock(r,-1);
    Range Arange,
          Brange;
    for(auto aio : dA.offsets)
        {
        //Compute bi, new block index of blk
        computeBlockInd(aio.block,Ais,Ablock);
        for(auto j : range(Ablock))
            Bblock.at(P.dest(j)) = Ablock[j];
        Arange.init(make_indexdim(Ais,Ablock));
        Brange.init(make_indexdim(Bis,Bblock));

        auto bblock = getBlock(dB,Bis,Bblock);
        auto bref = makeRef(bblock,&Brange);
        auto aref = makeTenRef(dA.data(),aio.offset,dA.size(),&Arange);

        bref += permute(aref,P);
        }
    }

template<typename T>
void
doTask(Order<IQIndex> const& O,
       QDense<T> & dB)
    {
    auto const dA = dB;
    auto Bis = O.is2();
    permuteQDense(O.perm(),dA,O.is1(),dB,Bis);
    }
template void doTask(Order<IQIndex> const&,QDense<Real> &);
template void doTask(Order<IQIndex> const&,QDense<Cplx> &);

} //namespace itensor

