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
#include "itensor/itdata/diag.h"
#include "itensor/itdata/qdiag.h"
#include "itensor/detail/gcounter.h"
#include "itensor/tensor/contract.h"
#include "itensor/itdata/qutil.h"
#include "itensor/util/iterate.h"

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
               IndexSet const& is,
               F const& callback)
    {
    auto r = order(is);
    auto block = IntArray(r,0);
    auto blockMinM = IntArray(r,0);
    auto blockSize = [&block,&is](long i)->long
        { 
        return (is[i]).blocksize0(block[i]); 
        };
    auto done = [&is,&block]()->bool
        {
        for(auto n : range(block)) 
            {
            if(block[n] >= is[n].nblock()) return true;
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
        auto ne = blockMinM[0]+blockSize(0);
        for(auto i : range(1,r))
            {
            ne = std::min(ne,blockMinM[i]+blockSize(i));
            }

        if(nb < ne) callback(nb,ne,block);

        //advance block indices and blockMinM
        for(auto i : range(r))
            {
            if(ne == blockMinM[i]+blockSize(i))
                {
                blockMinM[i] += blockSize(i);
                block[i] += 1;
                }
            }
        }
    }

template<typename T>
QN
doTask(CalcDiv const& C, QDiag<T> const& D)
    {
    if(order(C.is)==0) return QN();

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
        if(q != d) Error("Diagonal elements of QDiag ITensor would have inconsistent divergence");
        };
    loopDiagBlocks(D,C.is,checkBlock);
#endif

    return d;
    }
template QN doTask(CalcDiv const& C, QDiag<Real> const& D);
template QN doTask(CalcDiv const& C, QDiag<Cplx> const& D);

size_t
computeLength(IndexSet const& is)
    {
    if(order(is)==0) return 1ul;

    auto length = dim(is[0]);
    for(auto& I : is)
        {
        if(length != dim(I))
            Error("QDiag storage requires all IQIndices to be same size");
        }
    return length;
    }

template<typename T>
QDiag<T>::
QDiag(IndexSet const& is)
  : length(computeLength(is))
    {
    store.assign(length,0);
#ifdef DEBUG
    doTask(CalcDiv{is},*this);
#endif
    }
template QDiag<Real>::QDiag(IndexSet const& is);
template QDiag<Cplx>::QDiag(IndexSet const& is);

template<typename T>
QDiag<T>::
QDiag(IndexSet const& is, T val_)
  : val(val_),
    length(computeLength(is))
    {
#ifdef DEBUG
    doTask(CalcDiv{is},*this);
#endif
    }
template QDiag<Real>::QDiag(IndexSet const& is, Real val_);
template QDiag<Cplx>::QDiag(IndexSet const& is, Cplx val_);

template<typename T>
Cplx
doTask(GetElt& G, QDiag<T> const& D)
    {
    auto r = G.is.order();
#ifdef DEBUG
    if(G.is.order() != decltype(r)(G.inds.size())) 
        {
        printfln("is.order() = %d, ind.size() = %d",G.is.order(),G.inds.size());
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
template Cplx doTask(GetElt& G, QDiag<Real> const& D);
template Cplx doTask(GetElt& G, QDiag<Cplx> const& D);


template<typename T>
Cplx
doTask(SumEls, QDiag<T> const& d)
    {
    if(d.allSame())
        {
        return d.val*d.length;
        }

    T s = 0.;
    for(auto& el : d) s += el;
    return s;
    }
template Cplx doTask(SumEls, QDiag<Real> const& d);
template Cplx doTask(SumEls, QDiag<Cplx> const& d);

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
doTask(PrintIT& P, QDiag<T> const& d)
    {
    P.s << format("QDiag%s%s\n",typeName<T>(),d.allSame()?" (all same)":"");
    Real scalefac = 1.0;
    if(!P.x.isTooBigForReal()) scalefac = P.x.real0();
    else P.s << "(omitting too large scale factor)\n";

    auto ord = P.is.order();
    if(ord == 0) 
        {
        P.s << "  ";
        auto val = d.allSame() ? d.val : d.store.front();
        P.s << formatVal(scalefac*val) << "\n";
        return;
        }

    //callback function to pass to loopDiagBlocks:
    auto printBlock = [&P,&d,ord,scalefac]
        (size_t nb, 
         size_t ne,
         IntArray const& block)
        {
        //print indices for this block
        for(auto i : range(ord))
            {
            if(i > 0) P.s << ", ";
            P.s << P.is[i].blocksize0(block[i])
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
                for(auto ii : range1(ord))
                    {
                    P.s << (1+j);
                    if(ii < ord) P.s << ",";
                    }
                P.s << ") " << formatVal(val) << "\n";
                }
            }
        };

    loopDiagBlocks(d,P.is,printBlock);
    }
template void doTask(PrintIT& P, QDiag<Real> const& d);
template void doTask(PrintIT& P, QDiag<Cplx> const& d);

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
               IndexSet const& Dis,
               Labels const& DL,
               QDense<VT> const& T,
               IndexSet const& Tis,
               Labels const& TL,
               IndexSet const& Cis,
               Labels const& CL,
               ManageStore & m)
    {
    using VC = common_type<VT,VD>;
#ifdef DEBUG
    if(Dis.order() == 0) Error("QDiag order 0 case not handled");
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
            (DataRange<const VT> tblock, Block const& Tblockind,
             DataRange<const VD> dblock, Block const& Dblockind,
             DataRange<VC>       cblock, Block const& Cblockind)
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
            (DataRange<const VT> tblock, Block const& Tblockind,
             DataRange<const VD> dblock, Block const& Dblockind,
             DataRange<VC>       cblock, Block const& Cblockind)
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

template<typename T>
bool
isReplaceDelta(QDiag<T> const& d, IndexSet const& dis, Labels const& l)
    {
    if( (order(dis) == 2) && d.allSame() 
         && (d.val == 1.) && (dim(dis[0]) == dim(dis[1])) )
        {
        bool i1_is_contracted = (l[0] < 0);
        bool i2_is_contracted = (l[1] < 0);
        return (i1_is_contracted != i2_is_contracted);
        }
    return false;
    }

template<typename T1, typename T2>
void
doTask(Contract& C,
       QDiag<T1> const& d,
       QDense<T2> const& t,
       ManageStore& m)
    {
    Labels Lind,
           Rind,
           Nind;
    computeLabels(C.Lis,C.Lis.order(),C.Ris,C.Ris.order(),Lind,Rind);
    //TODO: add a case where there is a scaled delta function,
    //so the data also gets scaled
    if( isReplaceDelta(d,C.Lis,Lind) )
        {
        //println("doTask(Contract,Diag,Dense): isReplaceDelta = true");
        // We are contracting with a delta function that is replacing
        // a single index
        contractISReplaceIndex(C.Ris,Rind,C.Lis,Lind,C.Nis);

        // Output data is the dense storage
        m.makeNewData<QDense<T2>>(t.offsets,t.begin(),t.end());
        }
    else
        {
        bool sortIndices = false;
        contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
        blockDiagDense(d,C.Lis,Lind,
                       t,C.Ris,Rind,
                       C.Nis,Nind,m);
        }
    }
template void doTask(Contract& Con,QDiag<Real> const& A,QDense<Real> const& B,ManageStore& m);
template void doTask(Contract& Con,QDiag<Cplx> const& A,QDense<Real> const& B,ManageStore& m);
template void doTask(Contract& Con,QDiag<Real> const& A,QDense<Cplx> const& B,ManageStore& m);
template void doTask(Contract& Con,QDiag<Cplx> const& A,QDense<Cplx> const& B,ManageStore& m);

template<typename T1, typename T2>
void
doTask(Contract& C,
       QDense<T1> const& t,
       QDiag<T2> const& d,
       ManageStore& m)
    {
    Labels Lind,
           Rind,
           Nind;
    computeLabels(C.Lis,C.Lis.order(),C.Ris,C.Ris.order(),Lind,Rind);
    //TODO: add a case where there is a scaled delta function,
    //so the data also gets scaled
    if( isReplaceDelta(d,C.Ris,Rind) )
        {
        //println("doTask(Contract,QDense,QDiag): isReplaceDelta = true");
        // We are contracting with a delta function that is replacing
        // a single index
        contractISReplaceIndex(C.Lis,Lind,C.Ris,Rind,C.Nis);
        }
    else
        {
        bool sortInds = false;
        contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortInds);
        blockDiagDense(d,C.Ris,Rind,
                       t,C.Lis,Lind,
                       C.Nis,Nind,m);
        }
    }
template void doTask(Contract& Con,QDense<Real> const& A,QDiag<Real> const& B,ManageStore& m);
template void doTask(Contract& Con,QDense<Cplx> const& A,QDiag<Real> const& B,ManageStore& m);
template void doTask(Contract& Con,QDense<Real> const& A,QDiag<Cplx> const& B,ManageStore& m);
template void doTask(Contract& Con,QDense<Cplx> const& A,QDiag<Cplx> const& B,ManageStore& m);

template<typename T>
bool
doTask(IsDense,
       QDiag<T> const& d)
    {
    return false;
    }
template bool doTask(IsDense,QDiag<Real> const& d);
template bool doTask(IsDense,QDiag<Cplx> const& d);

template<typename V>
void
doTask(RemoveQNs & R, 
       QDiag<V> const& qd,
       ManageStore & m)
    {
    //auto *pd = qd.data();
    if(qd.allSame()) m.makeNewData<Diag<V>>(qd.length,qd.val);
    else             m.makeNewData<Diag<V>>(qd.begin(),qd.end());
    }
template void doTask(RemoveQNs &, QDiag<Real> const&, ManageStore &);
template void doTask(RemoveQNs &, QDiag<Cplx> const&, ManageStore &);

template<typename V>
void
doTask(ToDense & R,
       QDiag<V> const& d,
       ManageStore & m)
  {
  auto ddiv = doTask(CalcDiv{R.is},d);
  auto nd = m.makeNewData<QDense<V>>(R.is,ddiv);
  auto ninds = length(R.is);

  Range dense_range;
  for(auto const& io : nd->offsets)
    {
    // Check that the current block
    // is on the diagonal. If it is not,
    // no diagonal elements need to be set.
    auto diag_block = true;
    for(auto i : range(1,ninds))
      {
      if(io.block[i] != io.block[0])
        diag_block = false;
      }

    if(diag_block)
      {
      // Make a TenRef of the current block we want to assign to
      dense_range.init(make_indexdim(R.is,io.block));
      auto aref = makeTenRef(nd->data(),io.offset,nd->size(),&dense_range);

      long tot_stride = 0; //total strides for this block
      for(auto i : range(ninds))
        tot_stride += dense_range.stride(i);

      if( d.allSame() )
        {
        auto block_size = dense_range.extent(0); // Assume all dimensions of the block are
                                                 // the same (maybe this is not always true?)
        for(auto i : range(block_size))
          aref[i*tot_stride] = d.val;
        }
      else
        {
        // Get the diagonal elements of the current
        // block
        auto pD = getBlock(d,R.is,io.block);
        // Make a VecRef of the diagonal elements of the block to assign
        // to the diagonal of the new QDense storage block
        auto Dref = makeVecRef(pD.data(),pD.size());
        for(auto i : range(Dref.size()))
          aref[i*tot_stride] = Dref[i];
        }
      }
    }
  }
template void doTask(ToDense &, QDiag<Real> const&, ManageStore &);
template void doTask(ToDense &, QDiag<Cplx> const&, ManageStore &);

} //namespace itensor

