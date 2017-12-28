//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/print_macro.h"
#include "itensor/util/stdx.h"
#include "itensor/util/range.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/tensor/sliceten.h"
#include "itensor/iqtensor.h"
#include "itensor/itdata/qutil.h"

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
using std::move;



template<>
IQTensor::
ITensorT(Complex val) 
    :
    scale_(1.)
    { 
    if(val.imag()==0.)
        {
        store_ = newITData<ScalarReal>(val.real());
        }
    else
        {
        store_ = newITData<ScalarCplx>(val);
        }
    }

//IQTensor::
//IQTensor(const QN& q, vector<IQIndex>&& iqinds) 
//	: 
//    is_(move(iqinds)),
//    store_(newITData<IQTReal>(is_,q)),
//    div_(q),
//    scale_(1.)
//	{ }

//IQTensor::
//IQTensor(const QN& q,
//         IQIndexSet&& iset,
//         NewData nd,
//         LogNum scale)
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
//    Labels Lind,
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
    const Labels& block_ind;
    const Permutation& P;
    Real fac = 0;
    AddITensor(QN const& tdiv_,
               IQIndexSet const& iqis_,
               IndexSet const& is_,
               Labels const& block_ind_,
               Permutation const& P_,
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

const char*
typeNameOf(AddITensor const&) { return "AddITensor"; }

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
addIT(AddITensor & A, 
      QDense<T1> & d, 
      Dense<T2> const& t)
    {
    auto ddiv = doTask(CalcDiv{A.iqis},d);
    if(ddiv != A.tdiv) Error("IQTensor+=ITensor, ITensor has incompatible QN divergence");
    Range drange;
    drange.init(make_indexdim(A.iqis,A.block_ind));
    auto dblock = getBlock(d,A.iqis,A.block_ind);
    auto dref = makeRef(dblock,&drange);
    auto tref = makeTenRef(t.data(),t.size(),&A.is);
    transform(permute(tref,A.P),dref,Adder{A.fac});
    }

template<typename T1, typename T2>
void
doTask(AddITensor & A, 
       QDense<T1> const& d, 
       Dense<T2> const& t,
       ManageStore & m)
    {
    if(isReal(d) && isCplx(t))
        {
        auto *nd = m.makeNewData<QDenseCplx>(d.offsets,d.begin(),d.end());
        addIT(A,*nd,t);
        }
    else
        {
        auto *ncd = m.modifyData(d);
        addIT(A,*ncd,t);
        }
    }
template void doTask(AddITensor &, QDense<Real> const&, Dense<Real> const&, ManageStore &);
template void doTask(AddITensor &, QDense<Real> const&, Dense<Cplx> const&, ManageStore &);
template void doTask(AddITensor &, QDense<Cplx> const&, Dense<Real> const&, ManageStore &);
template void doTask(AddITensor &, QDense<Cplx> const&, Dense<Cplx> const&, ManageStore &);


IQTensor&
operator+=(IQTensor & T, 
           ITensor const& t)
    {
    if(!t) Error("IQTensor+=ITensor: r.h.s. ITensor is default constructed");
    if(!T.inds()) Error("Calling IQTensor+= on default constructed ITensor");
    auto rank = T.r();
#ifdef DEBUG
    if(t.r() != rank) Error("Mismatched number of indices in IQTensor+=ITensor");
#endif

    Permutation P(rank);
    Labels block_ind(rank);
    decltype(rank) nfound = 0;
    for(auto i : range(rank))
    for(auto I : range(rank))
        {
        auto j = findindex(T.inds()[I],t.inds()[i]);
        if(j > 0)
            {
            block_ind[I] = (j-1);
            P.setFromTo(i,I);
            ++nfound;
            break;
            }
        }

    if(nfound != rank)
        {
        Print(T.inds());
        Print(t.inds());
        Error("Adding IQTensor+=ITensor, index of ITensor not a subset of any IQIndex");
        }

    auto tdiv = calcDiv(T.inds(),block_ind);

    if(!T.store()) 
        {
        //allocate data to add this ITensor into
        if(!isComplex(t)) T.store() = newITData<QDenseReal>(T.inds(),tdiv);
        else              T.store() = newITData<QDenseCplx>(T.inds(),tdiv);
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


template<typename V>
ITensor
doTask(ToITensor & T, 
       QDense<V> const& d)
    {
    auto r = T.is.r();
    auto nd = Dense<V>(area(T.is),0);
    auto *pd = d.data();
    auto *pn = nd.data();
    vector<long> block(r,0);
    detail::GCounter C(r);
    for(auto& io : d.offsets)
        {
        computeBlockInd(io.block,T.is,block);
        for(long j = 0; j < r; ++j)
            {
            long start = 0;
            for(long b = 0; b < block[j]; ++b)
                start += T.is[j][b].m();
            C.setRange(j,start,start+T.is[j][block[j]].m()-1);
            }
        //TODO: need to make a Range/TensoRef iterator
        //to rewrite the following code more efficiently
        for(; C.notDone(); ++C)
            {
            pn[offset(T.is,C.i)] = pd[io.offset+C.ind];
            }
        }
    auto inds = IndexSetBuilder(r);
    for(decltype(r) j = 0; j < r; ++j) inds.setIndex(j,T.is[j]);
    return ITensor(inds.build(),std::move(nd),T.scale);
    }
template ITensor doTask(ToITensor & T, QDense<Real> const& d);
template ITensor doTask(ToITensor & T, QDense<Cplx> const& d);

template<typename V>
ITensor
doTask(ToITensor & T, QDiag<V> const& qd)
    {
    Diag<V> d;
    if(qd.allSame()) d = Diag<V>(qd.length,qd.val);
    else             d = Diag<V>(qd.begin(),qd.end());
    auto r = T.is.r();
    auto inds = IndexSetBuilder(r);
    for(auto j : range(r)) inds.setIndex(j,T.is[j]);
    return ITensor{inds.build(),std::move(d),T.scale};
    }
template ITensor doTask(ToITensor & T, QDiagReal const& d);
template ITensor doTask(ToITensor & T, QDiagCplx const& d);

template<typename V>
ITensor
doTask(ToITensor & T, 
       QMixed<V> const& d)
    {
    auto inds = IndexSetBuilder(rank(T.is));
    for(auto j : range(rank(T.is))) inds.setIndex(j,T.is[j]);
    return ITensor(inds.build(),Dense<V>(d.begin(),d.end()),T.scale);
    }
template ITensor doTask(ToITensor & T, QMixed<Real> const& d);
template ITensor doTask(ToITensor & T, QMixed<Cplx> const& d);

//template<typename D>
//ITensor
//doTask(ToITensor & T, ITDiag<D> const& d)
//    {
//    auto r = T.is.r();
//    IndexSet::storage_type inds(r);
//    for(long j = 0; j < r; ++j) inds[j].ext = T.is[j];
//    return ITensor(IndexSet{std::move(inds)},ITDiag<D>{d});
//    }
//template ITensor doTask(ToITensor & T, ITDiag<Real> const& d);
//template ITensor doTask(ToITensor & T, ITDiag<Cplx> const& d);


ITensor
toITensor(IQTensor const& T)
    {
    //Handle unallocated IQTensors
    if(!T.store()) 
        {
        if(rank(T)==0) return ITensor{};
        auto inds = IndexSetBuilder(rank(T));
        for(auto j : range(rank(T))) inds.setIndex(j,T.inds()[j]);
        return ITensor(inds.build());
        }
    //Main case for allocated IQTensors
    return doTask(ToITensor{T.inds(),T.scale()},T.store());
    }

QN
div(IQTensor const& T) 
    { 
    if(!T) Error("div(IQTensor) not defined for unallocated IQTensor");
    return doTask(CalcDiv{T.inds()},T.store());
    }

IQTensor
combiner(std::vector<IQIndex> cinds,
         Args const& args)
    {
    if(cinds.empty()) Error("No indices passed to combiner");
    auto cname = args.getString("IndexName","cmb");
    auto itype = getIndexType(args,"IndexType",Link);

    auto cdir = Out;
    if(args.defined("IndexDir"))
        {
        cdir = toArrow(args.getInt("IndexDir"));
        }
    else
        {
        //If not specified by user, make combined IQIndex
        //point Out unless all combined indices have In arrows.
        auto allin = true;
        for(auto& i : cinds) 
            if(i.dir() != In)
                {
                allin = false;
                break;
                }
        if(allin) cdir = In;
        }


    auto C = QCombiner{cinds};

    //Build the combined IQIndex,
    //saving information about
    //how we're doing it in C
    struct QNm { QN q; long m = 1; };
    auto qms = vector<QNm>{};

    //Loop over all possible QN sectors that
    //can be formed by the combined indices
    for(auto I : C.range())
        {
        QNm qm;
        //For this sector, figure out the total QN (qm.q)
        //and combined sector size (qm.m)
        for(auto j : range(cinds))
            {
            qm.q += cinds[j].qn(1+I[j]) * cinds[j].dir() * cdir;
            qm.m *= cinds[j].index(1+I[j]).m();
            }

        size_t n = 0;
        for(; n < qms.size(); ++n) if(qms[n].q==qm.q) break;

        if(n < qms.size())
            {
            C.setBlockRange(I,n,qms[n].m,qm.m);
            qms[n].m += qm.m;
            }
        else 
            {
            C.setBlockRange(I,n,0,qm.m);
            qms.push_back(qm);
            }
        }

    auto cstore = stdx::reserve_vector<IndexQN>(qms.size());
    for(auto n : range(qms)) 
        {
        cstore.emplace_back(Index{nameint("c",n),qms[n].m,itype},qms[n].q);
        }
    auto cind = IQIndex{cname,std::move(cstore),cdir};

    auto newind = IQIndexSetBuilder(1+cinds.size());
    newind.nextIndex(std::move(cind));
    for(auto& I : cinds) 
        {
        I.dag();
        newind.nextIndex(std::move(I));
        }

    return IQTensor{newind.build(),std::move(C)};
    }

struct IsQCombiner
    {
    template<typename D>
    bool 
    operator()(D const& d) { return false; }
    bool
    operator()(QCombiner const& d) { return true; }
    };

IQIndex
combinedIndex(IQTensor const& C)
    {
#ifdef DEBUG
    auto iscombiner = applyFunc(IsQCombiner{},C.store());
    if(not iscombiner)
        {
        throw ITError("Called combinedIndex on ITensor that is not a combiner");
        }
#endif
    return C.inds().front();
    }

IQIndex
findIQInd(IQTensor const& T, 
          Index    const& i)
    {
    for(auto& J : T.inds())
        if(hasindex(J,i)) return J;
    return IQIndex{};
    }

QN
qn(IQTensor const& T, Index const& i) 
    { 
    auto I = findIQInd(T,i);
    if(not I) Error("qn: no matching IQIndex found");
    return qn(I,i);
    }

Arrow
dir(IQTensor const& T, Index const& i) 
    { 
    auto I = findIQInd(T,i);
    if(not I) Error("dir: no matching IQIndex found");
    return I.dir(); 
    }


Arrow
dir(IQTensor const& T, 
    IQIndex  const& I)
	{
    for(const IQIndex& J : T.inds())
        if(I == J) return J.dir();
    Error("dir: IQIndex not found");
    return Out;
	}

bool
isEmpty(IQTensor const& T)
    {
    if(not T.store()) return true;
    return doTask(IsEmpty{},T.store());
    }

//template<typename T>
//void
//doTask(PrintIT<IQIndex>& P, const ITDiag<T>& d)
//    {
//    constexpr auto type = std::is_same<T,Real>::value ? "Real" : "Cplx";
//    P.printInfo(d,format("Diag %s%s",type,d.allSame()?", all same":""),
//                doTask(NormNoScale{},d));
//
//    auto r = P.is.r();
//
//    if(r == 0) 
//        {
//        P.s << "  ";
//        P.printVal(P.scalefac*(d.empty() ? d.val : d.store.front()));
//        return;
//        }
//
//    if(!P.print_data) return;
//
//    for(auto i : range(d.length))
//        {
//        auto val = P.scalefac*(d.allSame() ? d.val : d.store[i]);
//        if(std::norm(val) > Global::printScale())
//            {
//            P.s << "(";
//            for(decltype(r) j = 1; j < r; ++j)
//                {
//                P.s << (1+i) << ",";
//                }
//            P.s << (1+i) << ") ";
//            P.printVal(val);
//            }
//        }
//    }
//template void doTask(PrintIT<IQIndex>& P, const ITDiag<Real>& d);
//template void doTask(PrintIT<IQIndex>& P, const ITDiag<Cplx>& d);


std::ostream&
operator<<(std::ostream& s, const IQTensor& T)
    {
	s << "/--------------IQTensor--------------\n";
    if(T.store())
        {
        s << "r=" << rank(T);
        s << " div=";
        try { 
            s << div(T); 
            }
        catch(ITError const& e)
            {
            s << "(undef.)";
            }
        s << " log(scale)=" << T.scale().logNum() << "\n";
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
