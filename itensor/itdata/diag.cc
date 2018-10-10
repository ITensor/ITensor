#include "itensor/itdata/diag.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/range.h"

using std::vector;

namespace itensor {

const char*
typeNameOf(DiagReal const& d) { return "DiagReal"; }
const char*
typeNameOf(DiagCplx const& d) { return "DiagCplx"; }

template <typename T>
Cplx
doTask(GetElt<Index> const& g, Diag<T> const& d)
    {
    auto first_i = (g.inds.empty() ? 0 : g.inds.front());
    //Check if inds_ reference an
    //element on the diagonal, else zero
    for(auto i : g.inds) if(i != first_i) return 0;

    if(d.allSame()) return d.val;
    return d.store.at(first_i);
    }
template Cplx doTask(GetElt<Index> const&, DiagReal const&);
template Cplx doTask(GetElt<Index> const&, DiagCplx const&);

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

template<typename T1, typename T2>
void
contractDiagDense(Diag<T1>  const& d,
                  IndexSet  const& dis,
                  Labels     const& dind,
                  Dense<T2> const& t,
                  IndexSet  const& tis,
                  Labels     const& tind,
                  Labels     const& Nind,
                  IndexSet  const& Nis,
                  ManageStore    & m)
    {
    using T3 = common_type<T1,T2>;
    bool t_has_uncontracted = false;
    for(auto j : range(tind)) 
        if(tind[j] >= 0)
            {
            t_has_uncontracted = true;
            break;
            }

    auto Tref = makeTenRef(t.data(),t.size(),&tis);

    if(t_has_uncontracted)
        {
        auto nd = m.makeNewData<Dense<T3>>(area(Nis),0.);
        auto Nref = makeTenRef(nd->data(),nd->size(),&Nis);
        if(d.allSame())
            {
            auto dref = UnifVecWrapper<decltype(d.val)>(d.val,d.length);
            contractDiagPartial(dref,dind,
                                Tref,tind,
                                Nref,Nind);
            }
        else
            {
            auto dref = makeVecRefc(d.data(),d.size());
            contractDiagPartial(dref,dind,
                                Tref,tind,
                                Nref,Nind);
            }
        }
    else //all inds of t contracted with d
        {
        long d_ustride = 0; //total result-stride of uncontracted inds of d
        for(auto i : range(dind))
            {
            if(dind[i] >= 0) d_ustride += dis.stride(i);
            }

        size_t nsize = (d_ustride==0) ? 1 : d.length;
        using nstorage_type = typename Diag<T3>::storage_type;
        auto nstore = nstorage_type(nsize,0);
        auto Nref = makeVecRef(nstore.data(),nsize);

        if(d.allSame())
            {
            auto dref = UnifVecWrapper<decltype(d.val)>(d.val,d.length);
            contractDiagFull(dref,dind,
                             Tref,tind,
                             Nref,Nind);
            }
        else
            {
            auto dref = makeVecRef(d.data(),d.size());
            contractDiagFull(dref,dind,
                             Tref,tind,
                             Nref,Nind);
            }
        if(rank(Nis)==1)
            {
            m.makeNewData<Dense<T3>>(std::move(nstore));
            }
        else
            {
            if(nsize==1)
                m.makeNewData<Diag<T3>>(1,nstore.front());
            else
                m.makeNewData<Diag<T3>>(std::move(nstore));
            }
        }
    }

template<typename T1, typename T2>
void
doTask(Contract<Index> & C,
       Dense<T1>  const& t,
       Diag<T2>   const& d,
       ManageStore     & m)
    { 
    Labels Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    bool sortIndices = false;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
    contractDiagDense(d,C.Ris,Rind,t,C.Lis,Lind,Nind,C.Nis,m);
    }
template void doTask(Contract<Index>&, Dense<Real> const&, Diag<Real> const&, ManageStore&);
template void doTask(Contract<Index>&, Dense<Real> const&, Diag<Cplx> const&, ManageStore&);
template void doTask(Contract<Index>&, Dense<Cplx> const&, Diag<Real> const&, ManageStore&);
template void doTask(Contract<Index>&, Dense<Cplx> const&, Diag<Cplx> const&, ManageStore&);

template<typename T1, typename T2>
void
doTask(Contract<Index> & C,
       Diag<T1>   const& d,
       Dense<T2>  const& t,
       ManageStore     & m)
    {
    Labels Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    bool sortIndices = false;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
    contractDiagDense(d,C.Lis,Lind,t,C.Ris,Rind,Nind,C.Nis,m);
    }
template void doTask(Contract<Index>&, Diag<Real> const&, Dense<Real> const&, ManageStore&);
template void doTask(Contract<Index>&, Diag<Real> const&, Dense<Cplx> const&, ManageStore&);
template void doTask(Contract<Index>&, Diag<Cplx> const&, Dense<Real> const&, ManageStore&);
template void doTask(Contract<Index>&, Diag<Cplx> const&, Dense<Cplx> const&, ManageStore&);

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
add(PlusEQ<Index> const& P,
    Diag<T1>          & D1,
    Diag<T2>     const& D2)
    {
#ifdef DEBUG
    if(D1.length != D2.length) Error("Mismatched lengths in plusEq");
#endif
    if(D1.allSame() || D2.allSame()) Error("Diag plusEq allSame case not implemented");

    if(std::is_same<T1,T2>::value)
        {
        auto d1 = realData(D1);
        auto d2 = realData(D2);
        daxpy_wrapper(d1.size(),P.alpha(),d2.data(),1,d1.data(),1);
        }
    else
        {
        auto ref1 = makeVecRef(D1.data(),D1.size());
        auto ref2 = makeVecRef(D2.data(),D2.size());
        transform(ref2,ref1,Adder{P.alpha()});
        }
    }

template<typename T1, typename T2>
void
doTask(PlusEQ<Index> const& P,
       Diag<T1> const& D1,
       Diag<T2> const& D2,
       ManageStore & m)
    {
    if(isReal(D1) && isCplx(D2))
        {
        auto *ncD1 = m.makeNewData<DiagCplx>(D1.begin(),D1.end());
        add(P,*ncD1,D2);
        }
    else
        {
        auto *ncD1 = m.modifyData(D1);
        add(P,*ncD1,D2);
        }
    }
template void doTask(PlusEQ<Index> const&,Diag<Real> const&,Diag<Real> const&,ManageStore &);
template void doTask(PlusEQ<Index> const&,Diag<Real> const&,Diag<Cplx> const&,ManageStore &);
template void doTask(PlusEQ<Index> const&,Diag<Cplx> const&,Diag<Real> const&,ManageStore &);
template void doTask(PlusEQ<Index> const&,Diag<Cplx> const&,Diag<Cplx> const&,ManageStore &);

template<typename N, typename T>
void
doTask(Fill<N>  const& f, 
       Diag<T>  const& d, 
       ManageStore   & m)
    {
    m.makeNewData<Diag<N>>(d.length,f.x);
    }
template void doTask(Fill<Real> const& f, DiagReal const& d, ManageStore& m);
template void doTask(Fill<Real> const& f, DiagCplx const& d, ManageStore& m);
template void doTask(Fill<Cplx> const& f, DiagReal const& d, ManageStore& m);
template void doTask(Fill<Cplx> const& f, DiagCplx const& d, ManageStore& m);

template<typename T1, typename T2>
void
doMult(T1 fac, Diag<T2> & D)
    {
    D.val *= fac;
    for(auto& elt : D.store) elt *= fac;
    }

template<typename T1, typename T2>
void
doTask(Mult<T1> const& M, Diag<T2> & D)
    {
    doMult(M.x,D);
    }
template void doTask(Mult<Real> const&, Diag<Real> &);
template void doTask(Mult<Real> const&, Diag<Cplx> &);
template void doTask(Mult<Cplx> const&, Diag<Cplx> &);

void
doTask(Mult<Cplx> const& M, Diag<Real> const& D, ManageStore & m)
    {
    auto *nD = m.makeNewData<DiagCplx>(D);
    doMult(M.x,*nD);
    }

template<typename T>
Real
doTask(NormNoScale, Diag<T> const& D)
    {
    if(D.allSame()) return std::sqrt(std::norm(D.val))*std::sqrt(D.length);
    auto d = realData(D);
    return dnrm2_wrapper(d.size(),d.data());
    }
template Real doTask(NormNoScale, Diag<Real> const& d);
template Real doTask(NormNoScale, Diag<Cplx> const& d);

void
doTask(Conj, DiagReal const& d) { }

void
doTask(Conj, DiagCplx & d) 
    { 
    if(d.allSame()) 
        {
        d.val = std::conj(d.val);
        }
    else
        {
        for(auto& el : d.store) 
            el = std::conj(el);
        }
    }


void
doTask(TakeReal, DiagReal const& D) { }

void
doTask(TakeReal, DiagCplx const& D, ManageStore& m) 
    { 
    if(D.allSame()) 
        {
        m.makeNewData<DiagReal>(D.length,D.val.real());
        }
    else            
        {
        auto nD = m.makeNewData<DiagReal>(D.size());
        for(auto i : range(D.store)) nD->store[i] = D.store[i].real();
        }
    }

void
doTask(TakeImag, DiagReal & D)
    {
    D.val = 0.;
    for(auto& el : D.store) el = 0.;
    }

void
doTask(TakeImag, DiagCplx const& D, ManageStore& m) 
    { 
    if(D.allSame()) 
        {
        m.makeNewData<DiagReal>(D.length,D.val.imag());
        }
    else            
        {
        auto nD = m.makeNewData<DiagReal>(D.size());
        for(auto i : range(D.store)) nD->store[i] = D.store[i].imag();
        }
    }

template<typename T>
void
doTask(PrintIT<Index>& P, Diag<T> const& d)
    {
    auto type = std::is_same<T,Real>::value ? "Real" : "Cplx";
    P.printInfo(d,format("Diag %s%s",type,d.allSame()?", all same":""),
              doTask(NormNoScale{},d));

    auto r = P.is.r();

    if(r == 0) 
        {
        P.s << "  ";
        P.s << formatVal(P.scalefac*(d.empty() ? d.val : d.store.front())) << "\n";
        return;
        }

    if(!P.print_data) return;

    for(auto i : range(d.length))
        {
        auto val = P.scalefac*(d.allSame() ? d.val : d.store[i]);
        if(std::norm(val) >= Global::printScale())
            {
            P.s << "(";
            for(decltype(r) j = 1; j < r; ++j)
                {
                P.s << (1+i) << ",";
                }
            P.s << (1+i) << ") ";
            P.s << formatVal(val) << "\n";
            }
        }
    }
template void doTask(PrintIT<Index>& P, DiagReal const& d);
template void doTask(PrintIT<Index>& P, DiagCplx const& d);

template <class T>
Cplx
doTask(SumEls<Index> S, Diag<T> const& d) 
    { 
    if(d.allSame()) return Real(minM(S.is))*d.val;
    T sum = 0;
    for(const auto& elt : d.store)
        sum += elt;
    return sum;
    }
template Cplx doTask(SumEls<Index> S, DiagReal const& d);
template Cplx doTask(SumEls<Index> S, DiagCplx const& d);

} //namespace itensor
