#include "itensor/itdata/itdiag.h"
#include "itensor/tensor/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"

using std::vector;

namespace itensor {

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

class UnifVecWrapper
    {
    Real val_;
    size_t size_;
    public:
    UnifVecWrapper(Real v, size_t s) : val_(v), size_(s) { }

    size_t
    size() const { return size_; }

    Real
    operator()(size_t j) const { return val_; }
    };

void
contractDiagDense(DiagReal  const& d,
                  IndexSet  const& dis,
                  Label     const& dind,
                  DenseReal const& t,
                  IndexSet  const& tis,
                  Label     const& tind,
                  Label     const& Nind,
                  IndexSet  const& Nis,
                  ManageStore    & m)
    {
    bool t_has_uncontracted = false;
    for(auto j : index(tind)) 
        if(tind[j] >= 0)
            {
            t_has_uncontracted = true;
            break;
            }

    auto Tref = makeTenRef(t.data(),t.size(),&tis);

    if(t_has_uncontracted)
        {
        auto nd = m.makeNewData<DenseReal>(area(Nis),0.);
        auto Nref = makeTenRef(nd->data(),nd->size(),&Nis);
        if(d.allSame())
            {
            auto dref = UnifVecWrapper(d.val,d.length);
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
        for(auto i : index(dind))
            {
            if(dind[i] >= 0) d_ustride += dis.stride(i);
            }

        size_t nsize = (d_ustride==0) ? 1 : d.length;
        auto nstore = DiagReal::storage_type(nsize,0);
        auto Nref = makeVecRef(nstore.data(),nsize);

        if(d.allSame())
            {
            auto dref = UnifVecWrapper(d.val,d.length);
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
        if(nsize==1)
            m.makeNewData<DiagReal>(1,nstore.front());
        else
            m.makeNewData<DiagReal>(std::move(nstore));
        }
    }

void
doTask(Contract<Index> & C,
       DenseReal  const& t,
       DiagReal   const& d,
       ManageStore     & m)
    { 
    Label Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    bool sortIndices = false;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
    contractDiagDense(d,C.Ris,Rind,t,C.Lis,Lind,Nind,C.Nis,m);
    }
void
doTask(Contract<Index> & C,
       DiagReal   const& d,
       DenseReal  const& t,
       ManageStore     & m)
    {
    Label Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    bool sortIndices = false;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
    contractDiagDense(d,C.Lis,Lind,t,C.Ris,Rind,Nind,C.Nis,m);
    }

void
doTask(PlusEQ<Index> const& P,
       DiagReal           & D1,
       DiagReal      const& D2)
    {
#ifdef DEBUG
    if(D1.size() != D2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(D1.allSame() || D2.allSame()) Error("Diag plusEq allSame case not implemented");
    daxpy_wrapper(D1.size(),P.fac(),D2.data(),1,D1.data(),1);
    }

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
        for(auto i : index(D.store)) nD->store[i] = D.store[i].real();
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
        for(auto i : index(D.store)) nD->store[i] = D.store[i].imag();
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
