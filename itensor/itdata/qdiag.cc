//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/qdiag.h"
#include "itensor/detail/gcounter.h"
#include "itensor/tensor/contract.h"
#include "itensor/itdata/qutil.h"
#include "itensor/util/range.h"

using std::vector;
using std::move;

namespace itensor {

const char*
typeNameOf(QDiagReal const& d) { return "QDiagReal"; }
const char*
typeNameOf(QDiagCplx const& d) { return "QDiagCplx"; }

template<typename T, class F>
void
loopDiagBlocks(QDiag<T> const& D,
               IQIndexSet const& is,
               F const& callback)
    {
    auto r = rank(is);
    auto block = IntArray(r,0);
    auto blockMinM = IntArray(r,0);
    auto blockIndex = [&block,&is](long i)->Index 
        { 
        return (is[i])[block[i]]; 
        };
    auto done = [&is,&block]()->bool
        {
        for(auto n : range(block)) 
            {
            if(block[n] >= is[n].nindex()) return true;
            }
        return false;
        };

    while(not done())
        {
        //diag elems start at nb
        auto nb = blockMinM[0];
        for(auto i : range(1,r))
            {
            nb = std::max(nb,blockMinM[i]);
            }

        //diag elems stop before ne
        auto ne = blockMinM[0]+blockIndex(0).m();
        for(auto i : range(1,r))
            {
            ne = std::min(ne,blockMinM[i]+blockIndex(i).m());
            }

        if(nb < ne) callback(nb,ne,block);

        //advance block indices and blockMinM
        for(auto i : range(r))
            {
            if(ne == blockMinM[i]+blockIndex(i).m())
                {
                blockMinM[i] += blockIndex(i).m();
                block[i] += 1;
                }
            }
        }
    }

template<typename T>
QN
doTask(CalcDiv const& C, QDiag<T> const& D)
    {
    if(rank(C.is)==0) return QN();

    auto d = QN();
    for(auto n : range(C.is))
        {
        d += C.is[n].qn(1)*C.is[n].dir();
        }

#ifdef DEBUG
    auto checkBlock = [&d,&C]
        (size_t nb, 
         size_t ne,
         IntArray const& block)
        {
        auto q = QN();
        for(auto n : range(C.is))
            {
            q += C.is[n].qn(1+block[n])*C.is[n].dir();
            }
        if(q != d) Error("Diagonal elements of QDiag IQTensor would have inconsistent divergence");
        };
    loopDiagBlocks(D,C.is,checkBlock);
#endif

    return d;
    }
template QN doTask(CalcDiv const& C, QDiag<Real> const& D);
template QN doTask(CalcDiv const& C, QDiag<Cplx> const& D);

size_t
computeLength(IQIndexSet const& is)
    {
    if(rank(is)==0) return 1ul;

    auto length = is[0].m();
    for(auto& I : is)
        {
        if(length != I.m())
            Error("QDiag storage requires all IQIndices to be same size");
        }
    return length;
    }

template<typename T>
QDiag<T>::
QDiag(IQIndexSet const& is)
  : length(computeLength(is))
    {
    store.assign(length,0);
#ifdef DEBUG
    doTask(CalcDiv{is},*this);
#endif
    }
template QDiag<Real>::QDiag(IQIndexSet const& is);
template QDiag<Cplx>::QDiag(IQIndexSet const& is);

template<typename T>
QDiag<T>::
QDiag(IQIndexSet const& is, T val_)
  : val(val_),
    length(computeLength(is))
    {
#ifdef DEBUG
    doTask(CalcDiv{is},*this);
#endif
    }
template QDiag<Real>::QDiag(IQIndexSet const& is, Real val_);
template QDiag<Cplx>::QDiag(IQIndexSet const& is, Cplx val_);

template<typename T>
Cplx
doTask(GetElt<IQIndex>& G, QDiag<T> const& D)
    {
    auto r = G.is.r();
#ifdef DEBUG
    if(G.is.r() != decltype(r)(G.inds.size())) 
        {
        printfln("is.r() = %d, ind.size() = %d",G.is.r(),G.inds.size());
        Error("Wrong number of indices passed to .real or .cplx");
        }
#endif
    if(r == 0) return D.store.front();

    auto first_i = (G.inds.empty() ? 0 : G.inds.front());
    //Check if inds_ reference an
    //element on the diagonal, else zero
    for(auto i : G.inds) if(i != first_i) return 0;

    if(D.allSame()) return D.val;

    auto n = G.inds[0];
#ifdef DEBUG
    if(n > (decltype(n))D.size()) Error("index out of range in getElt(QDiag..)");
#endif
    for(auto i : range(1,r))
        {
        if(n != G.inds[i]) return 0.;
        }
    return D.store.at(n);
    }
template Cplx doTask(GetElt<IQIndex>& G, QDiag<Real> const& D);
template Cplx doTask(GetElt<IQIndex>& G, QDiag<Cplx> const& D);


template<typename T>
Cplx
doTask(SumEls<IQIndex>, QDiag<T> const& d)
    {
    if(d.allSame())
        {
        return d.val*d.length;
        }

    T s = 0.;
    for(auto& el : d) s += el;
    return s;
    }
template Cplx doTask(SumEls<IQIndex>, QDiag<Real> const& d);
template Cplx doTask(SumEls<IQIndex>, QDiag<Cplx> const& d);

template<typename T>
void
doTask(Mult<Real> const& M, QDiag<T>& D)
    {
    if(D.allSame())
        {
        D.val *= M.x;
        }
    else
        {
        auto d = realData(D);
        dscal_wrapper(d.size(),M.x,d.data());
        }
    }
template void doTask(Mult<Real> const&, QDiagReal&);
template void doTask(Mult<Real> const&, QDiagCplx&);

void
doTask(Mult<Cplx> const& M, QDiag<Cplx> & D)
    {
    if(D.allSame())
        {
        D.val *= M.x;
        }
    else
        {
        for(auto& el : D) el *= M.x;
        }
    }

void
doTask(Mult<Cplx> const& M, QDiag<Real> const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDiagCplx>();
    nd->length = d.length;
    nd->val = d.val;
    if(not d.allSame())
        {
        auto *nd = m.makeNewData<QDiagCplx>();
        nd->store = QDiagCplx::storage_type(d.begin(),d.end());
        }
    doTask(M,*nd);
    }

void
doTask(Conj, QDiagCplx & d)
    {
    if(d.allSame()) 
        {
        applyConj(d.val);
        }
    else
        {
        for(auto& el : d) applyConj(el);
        }
    }


template<typename T>
Real
doTask(NormNoScale, QDiag<T> const& D)
    { 
    if(D.allSame()) return std::sqrt(std::norm(D.val))*std::sqrt(D.length);
    auto d = realData(D);
    return dnrm2_wrapper(d.size(),d.data());
    }
template Real doTask(NormNoScale, QDiag<Real> const& D);
template Real doTask(NormNoScale, QDiag<Cplx> const& D);


template<typename T>
void
doTask(PrintIT<IQIndex>& P, QDiag<T> const& d)
    {
    P.s << format("QDiag%s%s\n",typeName<T>(),d.allSame()?" (all same)":"");
    Real scalefac = 1.0;
    if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    else P.s << "(omitting too large scale factor)\n";

    auto rank = P.is.r();
    if(rank == 0) 
        {
        P.s << "  ";
        auto val = d.allSame() ? d.val : d.store.front();
        P.s << formatVal(scalefac*val) << "\n";
        return;
        }

    //callback function to pass to loopDiagBlocks:
    auto printBlock = [&P,&d,rank,scalefac]
        (size_t nb, 
         size_t ne,
         IntArray const& block)
        {
        //print indices for this block
        for(auto i : range(rank))
            {
            if(i > 0) P.s << ", ";
            P.s << P.is[i][block[i]]
                << "<" << P.is[i].dir() << ">"
                << P.is[i].qn(1+block[i]);
            }
        P.s << "\n";

        //print diagonal elements in this block, if any
        for(auto j : range(nb,ne))
            {
            auto val = d.val;
            if(not d.allSame()) val = d.store.at(j);
            val *= scalefac;
            if(std::norm(val) >= Global::printScale())
                {
                P.s << "(";
                for(auto ii : range1(rank))
                    {
                    P.s << (1+j);
                    if(ii < rank) P.s << ",";
                    }
                P.s << ") " << formatVal(val) << "\n";
                }
            }
        };

    loopDiagBlocks(d,P.is,printBlock);
    }
template void doTask(PrintIT<IQIndex>& P, QDiag<Real> const& d);
template void doTask(PrintIT<IQIndex>& P, QDiag<Cplx> const& d);

template<typename T>
class UnifVecWrapper
    {
    T val_;
    size_t size_;
    public:
    UnifVecWrapper(T v, size_t s) : val_(v), size_(s) { }

    size_t
    size() const { return size_; }

    T
    operator()(size_t j) const { return val_; }
    };

template<typename VD, typename VT>
void
blockDiagDense(QDiag<VD> const& D,
               IQIndexSet const& Dis,
               Labels const& DL,
               QDense<VT> const& T,
               IQIndexSet const& Tis,
               Labels const& TL,
               IQIndexSet const& Cis,
               Labels const& CL,
               ManageStore & m)
    {
    using VC = common_type<VT,VD>;
#ifdef DEBUG
    if(Dis.r() == 0) Error("QDiag rank 0 case not handled");
#endif

    bool T_has_uncontracted = false;
    for(auto j : range(TL)) 
        if(TL[j] >= 0)
            {
            T_has_uncontracted = true;
            break;
            }

    auto Cdiv = doTask(CalcDiv{Dis},D)+doTask(CalcDiv{Tis},T);

    if(T_has_uncontracted)
        {
        auto *nd = m.makeNewData<QDense<VC>>(Cis,Cdiv);
        auto& C = *nd;

        auto do_contract =
            [&D,&Dis,&Tis,&Cis,&DL,&TL,&CL]
            (DataRange<const VT> tblock, IntArray const& Tblockind,
             DataRange<const VD> dblock, IntArray const& Dblockind,
             DataRange<VC>       cblock, IntArray const& Cblockind)
            {
            Range Trange,
                  Crange;
            Trange.init(make_indexdim(Tis,Tblockind));
            auto Tref = makeRef(tblock,&Trange);
            Crange.init(make_indexdim(Cis,Cblockind));
            auto Cref = makeRef(cblock,&Crange);

            long nb=-1,ne=-1;
            auto starts = IntArray{};
            std::tie(nb,ne,starts) = diagBlockBounds(Dis,Dblockind);
            assert(nb <= ne);
            auto Dsize = ne-nb;

            //println("In do_contract");
            //print("  ");Print(Tblockind);
            //print("  ");Print(Dblockind);
            //print("  ");Print(Cblockind);
            //print("  ");Print(Dsize);
            //print("  ");Print(starts);
            //print("  ");Print(dblock.size());
            //println();

            if(D.allSame())
                {
                auto dref = UnifVecWrapper<decltype(D.val)>(D.val,Dsize);
                contractDiagPartial(dref,DL,
                                    Tref,TL,
                                    Cref,CL,
                                    starts);
                }
            else
                {
                auto Dref = makeVecRef(dblock.data(),Dsize);
                contractDiagPartial(Dref,DL,
                                    Tref,TL,
                                    Cref,CL,
                                    starts);
                }
            };

        loopContractedBlocks(T,Tis,
                             D,Dis,
                             C,Cis,
                             do_contract);
        }
    else
        {
        auto *nd = m.makeNewData<QDiag<VC>>(Cis);
        auto& C = *nd;

        auto do_contract =
            [&D,&Dis,&Tis,&DL,&TL,&CL]
            (DataRange<const VT> tblock, IntArray const& Tblockind,
             DataRange<const VD> dblock, IntArray const& Dblockind,
             DataRange<VC>       cblock, IntArray const& Cblockind)
            {
            Range Trange;
            Trange.init(make_indexdim(Tis,Tblockind));
            auto Tref = makeRef(tblock,&Trange);

            long nb=-1,ne=-1;
            auto starts = IntArray{};
            std::tie(nb,ne,starts) = diagBlockBounds(Dis,Dblockind);
            assert(nb <= ne);
            auto Dsize = ne-nb;

            auto Cref = makeRef(cblock,VecRange(cblock.size()));

            if(D.allSame())
                {
                auto dref = UnifVecWrapper<decltype(D.val)>(D.val,Dsize);
                contractDiagFull(dref,DL,
                                 Tref,TL,
                                 Cref,CL,
                                 starts);
                }
            else
                {
                auto Dref = makeVecRef(dblock.data(),Dsize);
                contractDiagFull(Dref,DL,
                                 Tref,TL,
                                 Cref,CL,
                                 starts);
                }
            };

        loopContractedBlocks(T,Tis,
                             D,Dis,
                             C,Cis,
                             do_contract);
        }
    }

template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDiag<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m)
    {
    Labels AL,
           BL,
           CL;
    bool sortInds = false;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),AL,BL);
    contractIS(Con.Lis,AL,Con.Ris,BL,Con.Nis,CL,sortInds);
    blockDiagDense(A,Con.Lis,AL,
                   B,Con.Ris,BL,
                   Con.Nis,CL,m);
    }
template void doTask(Contract<IQIndex>& Con,QDiag<Real> const& A,QDense<Real> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDiag<Cplx> const& A,QDense<Real> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDiag<Real> const& A,QDense<Cplx> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDiag<Cplx> const& A,QDense<Cplx> const& B,ManageStore& m);

template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDense<VA> const& A,
       QDiag<VB> const& B,
       ManageStore& m)
    {
    Labels AL,
           BL,
           CL;
    bool sortInds = false;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),AL,BL);
    contractIS(Con.Lis,AL,Con.Ris,BL,Con.Nis,CL,sortInds);
    blockDiagDense(B,Con.Ris,BL,
                   A,Con.Lis,AL,
                   Con.Nis,CL,m);
    }
template void doTask(Contract<IQIndex>& Con,QDense<Real> const& A,QDiag<Real> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDense<Cplx> const& A,QDiag<Real> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDense<Real> const& A,QDiag<Cplx> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDense<Cplx> const& A,QDiag<Cplx> const& B,ManageStore& m);


} //namespace itensor

