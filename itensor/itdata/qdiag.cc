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
QDiag<T>::
QDiag(IQIndexSet const& is, 
      QN const& div)
    {
    auto totalsize = updateOffsets(is,div);
    store.assign(totalsize,0);
    }
template QDiag<Real>::QDiag(IQIndexSet const& is, QN const& div);
template QDiag<Cplx>::QDiag(IQIndexSet const& is, QN const& div);

template<typename T>
QDiag<T>::
QDiag(IQIndexSet const& is, 
      QN const& div,
      T val_)
  : val(val_)
    {
    updateOffsets(is,div);
    }
template QDiag<Real>::QDiag(IQIndexSet const& is, QN const& div, Real val_);
template QDiag<Cplx>::QDiag(IQIndexSet const& is, QN const& div, Cplx val_);

template<typename T>
long QDiag<T>::
updateOffsets(IQIndexSet const& is,
              QN const& div)
    {
    offsets.clear();

    if(is.r()==0)
        {
        offsets.push_back(make_blof(0,0));
        return 1;
        }

    auto C = detail::GCounter(is.r());
    for(auto j : range(is.r()))
        C.setRange(j,0,is[j].nindex()-1);

    long totalsize = 0;
    for(; C.notDone(); ++C)
        {
        QN blockqn;
        for(auto j : range(is.r()))
            {
            auto& J = is[j];
            blockqn += J.qn(1+C[j])*J.dir();
            }
        if(blockqn == div)
            {
            long indstr = 1, //accumulate Index strides
                 ind = 0,
                 minm = std::numeric_limits<long>::max();
            for(auto j : range(is.r()))
                {
                auto& J = is[j];
                auto i_j = C[j];
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
template long QDiag<Real>::updateOffsets(IQIndexSet const& is,QN const& div);
template long QDiag<Cplx>::updateOffsets(IQIndexSet const& is,QN const& div);

template<typename T>
Cplx
doTask(GetElt<IQIndex>& G, QDiag<T> const& D)
    {
    auto* pelt = getElt(D,G.is,G.inds);
    if(pelt) return *pelt;
    return 0.;
    }
template Cplx doTask(GetElt<IQIndex>& G, QDiag<Real> const& D);
template Cplx doTask(GetElt<IQIndex>& G, QDiag<Cplx> const& D);

template<typename T>
QN
doTask(CalcDiv const& C, QDiag<T> const& d)
    {
#ifdef DEBUG
    if(d.offsets.empty()) Error("Default constructed QDiag in doTask(CalcDiv,QDiag)");
#endif
    auto b = d.offsets.front().block;
    Labels block_ind(C.is.r());
    computeBlockInd(b,C.is,block_ind);
    return calcDiv(C.is,block_ind);
    }
template QN doTask(CalcDiv const& C, QDiag<Real> const& d);
template QN doTask(CalcDiv const& C, QDiag<Cplx> const& d);

template<typename T>
Cplx
doTask(SumEls<IQIndex>, QDiag<T> const& d)
    {
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
    auto d = realData(D);
    dscal_wrapper(d.size(),M.x,d.data());
    }
template void doTask(Mult<Real> const&, QDiagReal&);
template void doTask(Mult<Real> const&, QDiagCplx&);

void
doTask(Mult<Cplx> const& M, QDiag<Cplx> & d)
    {
    for(auto& el : d) el *= M.x;
    }

void
doTask(Mult<Cplx> const& M, QDiag<Real> const& d, ManageStore & m)
    {
    auto *nd = m.makeNewData<QDiagCplx>(d.offsets,d.begin(),d.end());
    doTask(M,*nd);
    }

void
doTask(Conj, QDiagCplx & d)
    {
    for(auto& el : d) applyConj(el);
    }


template<typename T>
Real
doTask(NormNoScale, QDiag<T> const& D)
    { 
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
        
    Labels block(rank,0);
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
            for(decltype(Ddim.size()) j = 0; j < Ddim.size(); ++j)
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

