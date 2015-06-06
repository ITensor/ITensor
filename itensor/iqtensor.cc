//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/count.h"
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/iqtensor.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::find;
using std::pair;
using std::make_pair;
using std::make_shared;
using std::shared_ptr;
using std::move;




template<> IQTensor::
ITensorT(const IQIndex& i1) :
    is_(i1),
    scale_(1.)
    { }

template<> IQTensor::
ITensorT(const IQIndex& i1,
         const IQIndex& i2) :
    is_(i1,i2),
    scale_(1.)
    { }


template<>
IQTensor::
ITensorT(Complex val) 
    :
    scale_(1.)
    { 
    if(val.imag()==0)
        store_ = make_shared<ITDiag<Real>>(1,val.real());
    else
        store_ = make_shared<ITDiag<Complex>>(1,val);
    }

//IQTensor::
//IQTensor(const QN& q, vector<IQIndex>&& iqinds) 
//	: 
//    is_(move(iqinds)),
//    store_(make_shared<IQTData>(is_,q)),
//    div_(q),
//    scale_(1.)
//	{ }

//IQTensor::
//IQTensor(const QN& q,
//         IQIndexSet&& iset,
//         NewData nd,
//         LogNumber scale)
//    :
//    is_(move(iset)),
//    store_(move(nd)),
//    div_(q),
//    scale_(scale)
//    {
//    }



//IQTensor&
//operator*=(IQTensor& A, const IQTensor& B)
//    {
//    if(!A || !B)
//        Error("Default constructed IQTensor in product");
//
//    if(&A == &B)
//        {
//        A = IQTensor(sqr(norm(A)));
//        return A;
//        }
//
//    auto& Lis = A.inds();
//    auto& Ris = B.inds();
//
//    auto checkDirs = 
//    [&Lis,&Ris](const IQIndex& li, const IQIndex& ri)
//        {
//        if(li.dir() == ri.dir())
//            {
//            println("-------------------------");
//            println("Left indices = \n",Lis);
//            println("-------------------------");
//            println("Right indices = \n",Ris);
//            println("-------------------------");
//            println("IQIndex from left = \n",li);
//            println("IQIndex from right = \n",ri);
//            Error("Incompatible arrow directions in IQTensor contraction");
//            }
//        };
//    Label Lind,
//          Rind;
//    computeLabels(Lis,Lis.r(),Ris,Ris.r(),Lind,Rind,checkDirs);
//
//    auto nstore = A.store();
//    auto C = doTask(Contract<IQIndex>{Lis,Lind,Ris,Rind},nstore,B.store());
//
//    auto nscale = A.scale()*B.scale();
//    if(!std::isnan(C.scalefac)) nscale *= C.scalefac;
//
//#ifdef DEBUG
//    //Check for duplicate indices
//    detail::check(C.Nis);
//#endif
//
//    A = IQTensor(C.Nis,std::move(nstore),nscale);
//
//    return A;
//    }

struct AddITensor
    {
    const QN& tdiv;
    const IQIndexSet& iqis;
    const IndexSet& is;
    const Label& block_ind;
    const Permutation& P;
    Real fac = 0;
    AddITensor(const QN& tdiv_,
               const IQIndexSet& iqis_,
               const IndexSet& is_,
               const Label& block_ind_,
               const Permutation& P_,
               Real scalefac_)
        :
        tdiv(tdiv_),
        iqis(iqis_),
        is(is_),
        block_ind(block_ind_),
        P(P_),
        fac(scalefac_)
        { }
    };

void
doTask(AddITensor& A, IQTData& d, const ITReal& t)
    {
    auto ddiv = calcDiv(A.iqis,d);
    if(ddiv != A.tdiv) Error("IQTensor+=ITensor, ITensor has incompatible QN flux/divergence");
    Range drange;
    drange.init(make_indexdim(A.iqis,A.block_ind));
    auto* dblock = d.getBlock(A.iqis,A.block_ind);

    auto dref = makeTensorRef(dblock,drange);
    auto tref = makeTensorRef(t.data(),A.is);
    auto add = [f=A.fac](Real& r1, Real r2) { r1 += f*r2; };
    permute(tref,A.P,dref,add);
    }


IQTensor&
operator+=(IQTensor& T, const ITensor& t)
    {
    if(!t) Error("IQTensor+=ITensor: r.h.s. ITensor is default constructed");
    if(!T.inds()) Error("Calling IQTensor+= on default constructed ITensor");
    auto rank = T.r();
#ifdef DEBUG
    if(t.r() != rank) Error("Mismatched number of indices in IQTensor+=ITensor");
#endif

    Permutation P(rank);
    Label block_ind(rank);
    for(auto i : count(rank))
    for(auto I : count(rank))
        {
        auto j = findindex(T.inds()[I],t.inds()[i]);
        if(j > 0)
            {
            block_ind[I] = (j-1);
            P.setFromTo(i,I);
            break;
            }
        }

    auto tdiv = calcDiv(T.inds(),block_ind);

    if(!T.store()) 
        {
        //allocate data to add this ITensor into
        if(!isComplex(t)) T.store() = make_shared<IQTData>(T.inds(),tdiv);
        else              Error("Initializing complex IQTensor in +=ITensor not yet implemented");
        }

    Real scalefac = 1;
    if(T.scale().magnitudeLessThan(t.scale())) T.scaleTo(t.scale()); 
    else                                       scalefac = (t.scale()/T.scale()).real();

    doTask(AddITensor{tdiv,T.inds(),t.inds(),block_ind,P,scalefac},T.store(),t.store());

    return T;
    }

template<>
IQTensor& IQTensor::
dag()
    {
    is_.dag();
    doTask(Conj{},store_);
    return *this;
    }


struct ToITensor
    {
    const IQIndexSet& is;
    const LogNumber& scale;

    ToITensor(const IQIndexSet& is_,
              const LogNumber& scale_)
        :
        is(is_),
        scale(scale_)
        { }
    };

ITensor
doTask(ToITensor& T, const IQTData& d)
    {
    auto r = T.is.r();
    auto nd = ITReal(area(T.is),0);
    auto *pd = d.data.data();
    auto *pn = nd.data();
    vector<long> block(r,0);
    detail::GCounter C(r);
    for(const auto& io : d.offsets)
        {
        inverseBlockInd(io.block,T.is,block);
        for(long j = 0; j < r; ++j)
            {
            long start = 0;
            for(long b = 0; b < block[j]; ++b)
                start += T.is[j][b].m();
            C.setInd(j,start,start+T.is[j][block[j]].m()-1);
            }
        //TODO: need to make a Range/TensoRef iterator
        //to rewrite the following code more efficiently
        for(; C.notDone(); ++C)
            {
            pn[ind(T.is,C.i)] = pd[io.offset+C.ind];
            }
        }
    IndexSet::storage_type inds(r);
    for(long j = 0; j < r; ++j) inds[j].ext = T.is[j];
    return ITensor(IndexSet{std::move(inds)},std::move(nd),T.scale);
    }

template<typename D>
ITensor
doTask(ToITensor& T, const ITDiag<D>& d)
    {
    auto r = T.is.r();
    IndexSet::storage_type inds(r);
    for(long j = 0; j < r; ++j) inds[j].ext = T.is[j];
    return ITensor(IndexSet{std::move(inds)},ITDiag<D>{d});
    }
template ITensor doTask(ToITensor& T, const ITDiag<Real>& d);
template ITensor doTask(ToITensor& T, const ITDiag<Cplx>& d);


ITensor
toITensor(const IQTensor& T)
    {
    //Handle unallocated IQTensors
    if(!T.store()) 
        {
        if(T.r()==0) return ITensor{};
        IndexSet::storage_type inds(T.r());
        for(long j = 0; j < T.r(); ++j) inds[j].ext = T.inds()[j];
        return ITensor(IndexSet{std::move(inds)});
        }
    //Main case for allocated IQTensors
    return doTask<ITensor>(ToITensor{T.inds(),T.scale()},T.store());
    }

struct CalcDiv 
    { 
    const IQIndexSet& is;
    CalcDiv(const IQIndexSet& is_) : is(is_) { }
    };

QN
doTask(const CalcDiv& C, const IQTData& d)
    {
    return calcDiv(C.is,d);
    }

QN
doTask(const CalcDiv& C, const ITCombiner& d)
    {
    return QN{};
    }

template<typename T>
QN
doTask(const CalcDiv& C, const ITDiag<T>& d)
    {
    return QN{};
    }
template QN doTask(const CalcDiv& C, const ITDiag<Real>& d);
template QN doTask(const CalcDiv& C, const ITDiag<Cplx>& d);

QN
div(const IQTensor& T) 
    { 
    if(!T) Error("div(IQTensor) not defined for unallocated IQTensor");
    return doTask<QN>(CalcDiv{T.inds()},T.store());
    }

IQTensor
combiner(std::vector<IQIndex> inds,
         const Args& args)
    {
    using QNm = pair<QN,long>;
    if(inds.empty()) Error("No indices passed to combiner");
    auto cname = args.getString("IndexName","cmb");
    auto itype = getIndexType(args,"IndexType",inds.front().type());
    auto r = inds.size();
    auto dir = inds.front().dir();

    //create combined index
    detail::GCounter C(r);
    long nsubblocks = 1;
    for(size_t j = 0; j < r; ++j)
        {
        C.setInd(j,1,inds[j].nindex());
        nsubblocks *= inds[j].nindex();
        }
    vector<QNm> qns(nsubblocks);
    for(auto& qm : qns)
        {
        assert(C.notDone());
        qm.second = 1;
        for(size_t j = 0; j < r; ++j)
            {
            qm.first += inds[j].qn(C.i[j])*inds[j].dir();
            qm.second *= inds[j].index(C.i[j]).m();
            }
        ++C;
        }
    //Sort qns by quantum number; there will be duplicates in general:
    auto qnless = [](const QNm& qm1, const QNm& qm2) { return qm1.first < qm2.first; };
    std::sort(qns.begin(),qns.end(),qnless);

    //Count number of sectors
    auto currqn = qns.front().first;
    long sectors = 1;
    for(const auto& qm : qns)
        {
        if(qm.first != currqn)
            {
            ++sectors;
            currqn = qm.first;
            }
        }
    IQIndex::storage indqn(sectors);
    currqn = qns.front().first;
    long currm = 0,
         count = 0;
    for(const auto& qm : qns)
        {
        if(qm.first != currqn)
            {
            indqn[count] = IndexQN(Index(nameint("C",count),currm,itype),currqn);
            currqn = qm.first;
            ++count;
            currm = qm.second;
            }
        else
            {
            currm += qm.second;
            }
        }
    //Handle last element of indqn:
    indqn[count] = IndexQN(Index(nameint("C",count),currm,itype),currqn);
    //Finally make new combined IQIndex and put at front of inds
    vector<IQIndex> newinds(inds.size()+1);
    auto it = newinds.begin();
    *it = IQIndex(cname,std::move(indqn),dir);
    for(const auto& I : inds)
        {
        ++it;
        *it = dag(I);
        }

    return IQTensor(IQIndexSet{std::move(newinds)},ITCombiner());
    }

IQIndex
findIQInd(const IQTensor& T, const Index& i)
    {
    for(const IQIndex& J : T.inds())
        if(hasindex(J,i)) return J;
    Print(T.inds());
    Print(i);
    throw ITError("Index i not found in any of T's IQIndices");
    return IQIndex{};
    }


Arrow
dir(const IQTensor& T, const IQIndex& I)
	{
    for(const IQIndex& J : T.inds())
        if(I == J) return J.dir();
    Error("dir: IQIndex not found");
    return Out;
	}

bool
isZero(const IQTensor& T, const Args& args)
    {
    Error("Not implemented");
    //if(T.empty()) return true;
    ////done with all fast checks
    //if(args.getBool("Fast",false)) return false;
    //for(const ITensor& t : T.blocks())
    //    {
    //    if(!isZero(t)) return false;
    //    }
    return true;
    }

template<typename T>
void
doTask(PrintIT<IQIndex>& P, const ITDiag<T>& d)
    {
    constexpr auto type = std::is_same<T,Real>::value ? "Real" : "Cplx";
    P.printInfo(d,format("Diag %s%s",type,d.allSame()?", all same":""),
                doTask(NormNoScale{},d));

    auto r = P.is.r();

    if(r == 0) 
        {
        P.s << "  ";
        P.printVal(P.scalefac*(d.empty() ? d.val : d.store.front()));
        return;
        }

    if(!P.print_data) return;

    for(auto i : count(d.length))
        {
        auto val = P.scalefac*(d.allSame() ? d.val : d.store[i]);
        if(std::norm(val) > Global::printScale())
            {
            P.s << "(";
            for(decltype(r) j = 1; j < r; ++j)
                {
                P.s << (1+i) << ",";
                }
            P.s << (1+i) << ") ";
            P.printVal(val);
            }
        }
    }
template void doTask(PrintIT<IQIndex>& P, const ITDiag<Real>& d);
template void doTask(PrintIT<IQIndex>& P, const ITDiag<Cplx>& d);


std::ostream&
operator<<(std::ostream& s, const IQTensor& T)
    {
	s << "/--------------IQTensor--------------\n";
    if(T.store())
        {
        s << "r=" << T.r() << " div=" << div(T) << " log(scale)=" << T.scale().logNum() << "\n";
        s << T.inds();
        //Checking whether std::ios::floatfield is set enables 
        //printing the contents of an ITensor when using the printf
        //format string %f (or another float-related format string)
        bool ff_set = (std::ios::floatfield & s.flags()) != 0;

        if(ff_set || Global::printdat())
            {
            s << "\n|-- Data -------\n";
            doTask(PrintIT<IQIndex>{s,T.scale(),T.inds(),true},T.store());
            s << "\n";
            }
        }
    else
        {
        s << "r=" << T.r() << " log(scale)=" << T.scale().logNum() << "\n";
        s << T.inds();
        s << "{Zero / Not yet allocated}\n";
        }
	s << "\\------------------------------------\n\n";
    return s;
    }


} //namespace itensor
