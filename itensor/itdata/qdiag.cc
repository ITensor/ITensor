//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/qdiag.h"
#include "itensor/detail/gcounter.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/range.h"

using std::vector;
using std::move;

namespace itensor {

template<typename T>
QN
doTask(CalcDiv const& C, QDiag<T> const& d)
    {
    if(rank(C.is)==0) return QN();

    auto d = QN();
    for(auto n : range(C.is))
        {
        d += C.is[n].qn(1)*C.is[n].dir();
        }

#ifdef DEBUG
    for(auto s : range(C.is[0].nindex()))
        {
        auto q = QN();
        for(auto n : range(C.is))
            {
            q += C.is[n].qn(1)*C.is[n].dir();
            }
        if(q != d) Error("Diagonal elements of QDiag IQTensor would have inconsistent flux");
        }
#endif

    return d;
    }
template QN doTask(CalcDiv const& C, QDiag<Real> const& d);
template QN doTask(CalcDiv const& C, QDiag<Cplx> const& d);

size_t
computeLength(IQIndexSet const& is)
    {
    if(rank(is)==0) return 1ul;

    auto length = is[0].m();
    for(auto& I : is)
        {
        if(length != is.m())
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
    if(D.allSame()) return D.val;

    auto r = G.is.r();
#ifdef DEBUG
    if(G.is.r() != decltype(r)(G.ind.size())) 
        {
        printfln("is.r() = %d, ind.size() = %d",G.is.r(),G.ind.size());
        Error("Wrong number of indices passed to .real or .cplx");
        }
#endif
    if(r == 0) return *(D.data());
    size_t n = G.ind[0];
#ifdef DEBUG
    if(n > D.size()) Error("index out of range in getElt(QDiag..)");
#endif
    for(auto i : range(1,r))
        {
        if(G.is[i]!=n) return 0.;
        }
    return *(D.data()+(n-1));
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
doTask(Mult<Cplx> const& M, QDiag<Cplx> & d)
    {
    if(D.allSame())
        {
        D.val *= M.x;
        }
    else
        {
        for(auto& el : d) el *= M.x;
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
    P.s << format("QDiag%s%s",typeName<T>(),d.allSame()?" (all same)":"");
    P.s << format(" {%d blocks; data size %d}\n",d.offsets.size(),d.size());
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
        
    auto block = IntArray(rank,0);
    auto blockIndex = [&block,&P](long i)->Index { return (P.is[i])[block[i]]; };

    Range brange;
    for(auto& io : d.offsets)
        {
        computeBlockInd(io.block,P.is,block);
        auto blockm = blockIndex(0).m();
        for(auto i : range(1,rank)) blockm = std::min(blockm,blockIndex(i).m());

        bool indices_printed = false;
        auto os = io.offset;
        for(auto n : range(blockm))
            {
            auto val = d.allSame() ? d.val : d.store[os++];
            val *= scalefac;
            if(std::norm(val) >= Global::printScale())
                {
                if(!indices_printed)
                    {
                    indices_printed = true;
                    //Print Indices of this block
                    for(auto i : range(rank))
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
                P.s << formatVal(val) << "\n";
                }
            }
        }
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
               Labels const& Dind,
               QDense<VT> const& T,
               IQIndexSet const& Tis,
               Labels const& Tind,
               IQIndexSet const& Cis,
               Labels const& Cind,
               ManageStore & m)
    {
    using VC = common_type<VT,VD>;
#ifdef DEBUG
    if(Dis.r() == 0) Error("QDiag rank 0 case not handled");
#endif

    bool T_has_uncontracted = false;
    for(auto j : range(Tind)) 
        if(Tind[j] >= 0)
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
            [&D,&Dis,&Tis,&Cis,&Dind,&Tind,&Cind]
            (DataRange<const VD> dblock, Labels const& Dblockind,
             DataRange<const VT> tblock, Labels const& Tblockind,
             DataRange<VC>       cblock, Labels const& Cblockind)
            {
            Range Trange,
                  Crange;
            Trange.init(make_indexdim(Tis,Tblockind));
            auto Tref = makeRef(tblock,&Trange);
            Crange.init(make_indexdim(Cis,Cblockind));
            auto Cref = makeRef(cblock,&Crange);

            auto Ddim = make_indexdim(Dis,Dblockind);
            auto Dminm = std::numeric_limits<size_t>::max();
            for(auto j : range(Ddim))
                {
                Dminm = std::min(Dminm,Ddim[j]);
                }

            if(D.allSame())
                {
                auto dref = UnifVecWrapper<decltype(D.val)>(D.val,Dminm);
                contractDiagPartial(dref,Dind,
                                    Tref,Tind,
                                    Cref,Cind);
                }
            else
                {
                auto Dref = makeVecRef(dblock.data(),Dminm);
                contractDiagPartial(Dref,Dind,
                                    Tref,Tind,
                                    Cref,Cind);
                }
            };

        loopContractedBlocks(D,Dis,
                             T,Tis,
                             C,Cis,
                             do_contract);
        }
    else
        {
        Error("Fully contracted QDiag*QDense not yet implemented");
        //auto nd = m.makeNewData<QDiag>(Cis,Cdiv);
        //auto& C = *nd;

        //loopContractedBlocks(D,Dis,
        //                     T,Tis,
        //                     C,Cis,
        //                     do_contract);
        }
    }

template<typename VA, typename VB>
void
doTask(Contract<IQIndex>& Con,
       QDiag<VA> const& A,
       QDense<VB> const& B,
       ManageStore& m)
    {
    Labels Aind,
           Bind,
           Cind;
    bool sortInds = false;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),Aind,Bind);
    contractIS(Con.Lis,Aind,Con.Ris,Bind,Con.Nis,Cind,sortInds);
    blockDiagDense(A,Con.Lis,Aind,
                   B,Con.Ris,Bind,
                   Con.Nis,Cind,m);
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
    Labels Aind,
           Bind,
           Cind;
    bool sortInds = false;
    computeLabels(Con.Lis,Con.Lis.r(),Con.Ris,Con.Ris.r(),Aind,Bind);
    contractIS(Con.Lis,Aind,Con.Ris,Bind,Con.Nis,Cind,sortInds);
    blockDiagDense(B,Con.Ris,Bind,
                   A,Con.Lis,Aind,
                   Con.Nis,Cind,m);
    }
template void doTask(Contract<IQIndex>& Con,QDense<Real> const& A,QDiag<Real> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDense<Cplx> const& A,QDiag<Real> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDense<Real> const& A,QDiag<Cplx> const& B,ManageStore& m);
template void doTask(Contract<IQIndex>& Con,QDense<Cplx> const& A,QDiag<Cplx> const& B,ManageStore& m);

} //namespace itensor

