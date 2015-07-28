#include "itensor/itdata/itdiag.h"
#include "itensor/itdata/itdata.h"
#include "itensor/matrix/lapack_wrap.h"
#include "itensor/tensor/contract.h"
#include "itensor/util/count.h"

using std::vector;

namespace itensor {

template <typename T>
Cplx
doTask(const GetElt<Index>& g, const ITDiag<T>& d)
    {
    auto first_i = (g.inds.empty() ? 0 : g.inds.front());
    //Check if inds_ reference an
    //element on the diagonal, else zero
    for(auto i : g.inds) if(i != first_i) return 0;

    if(d.allSame()) return d.val;
    return d.store.at(first_i);
    }
template Cplx doTask(const GetElt<Index>& g, const ITDiag<Real>& d);
template Cplx doTask(const GetElt<Index>& g, const ITDiag<Cplx>& d);

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
diagDense(ITDiag<Real> const& d,
          IndexSet const& dis,
          Label const& dind,
          ITReal const& t,
          IndexSet const& tis,
          Label const& tind,
          Label const& Nind,
          IndexSet & Nis,
          ManageStore & m)
    {
    bool t_has_uncontracted = false;
    for(auto j : index(tind)) 
        if(tind[j] >= 0)
            {
            t_has_uncontracted = true;
            break;
            }

    auto Tref = makeTensorRef(t.data(),tis);

    if(t_has_uncontracted)
        {
        auto nd = m.makeNewData<ITReal>(area(Nis),0.);
        auto Nref = makeTensorRef(nd->data(),Nis);
        if(d.allSame())
            {
            auto dref = UnifVecWrapper(d.val,d.length);
            contractDiagPartial(dref,dind,
                                Tref,tind,
                                Nref,Nind);
            }
        else
            {
            auto dref = VecRefc(d.data(),d.size());
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
        ITDiag<Real>::storage_type nstore(nsize,0);
        auto Nref = VecRef(nstore.data(),nsize);

        if(d.allSame())
            {
            auto dref = UnifVecWrapper(d.val,d.length);
            contractDiagFull(dref,dind,
                             Tref,tind,
                             Nref,Nind);
            }
        else
            {
            auto dref = VecRefc(d.data(),d.size());
            contractDiagFull(dref,dind,
                             Tref,tind,
                             Nref,Nind);
            }
        if(nsize==1)
            m.makeNewData<ITDiag<Real>>(1,nstore.front());
        else
            m.makeNewData<ITDiag<Real>>(std::move(nstore));
        }
    }

void
doTask(Contract<Index> & C,
       ITReal const& t,
       ITDiag<Real> const& d,
       ManageStore & m)
    { 
    Label Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    bool sortIndices = false;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
    diagDense(d,C.Ris,Rind,t,C.Lis,Lind,Nind,C.Nis,m);
    }
void
doTask(Contract<Index> & C,
       ITDiag<Real> const& d,
       ITReal const& t,
       ManageStore & m)
    {
    Label Lind,
          Rind,
          Nind;
    computeLabels(C.Lis,C.Lis.r(),C.Ris,C.Ris.r(),Lind,Rind);
    bool sortIndices = false;
    contractIS(C.Lis,Lind,C.Ris,Rind,C.Nis,Nind,sortIndices);
    diagDense(d,C.Lis,Lind,t,C.Ris,Rind,Nind,C.Nis,m);
    }

void
doTask(const PlusEQ<Index>& P,
       ITDiag<Real>& a1,
       const ITDiag<Real>& a2)
    {
#ifdef DEBUG
    if(a1.size() != a2.size()) Error("Mismatched sizes in plusEq");
#endif
    if(a1.allSame() || a2.allSame()) Error("ITDiag plusEq allSame case not implemented");
    daxpy_wrapper(a1.size(),P.fac,a2.data(),1,a1.data(),1);
    }

template<typename T>
void
doTask(const FillReal& f, const ITDiag<T>& d, ManageStore& m)
    {
    m.makeNewData<ITDiag<Real>>(d.length,f.r);
    }
template void doTask(const FillReal& f, const ITDiag<Real>& d, ManageStore& m);
template void doTask(const FillReal& f, const ITDiag<Cplx>& d, ManageStore& m);

template<typename T>
void
doTask(const FillCplx& f, const ITDiag<T>& d, ManageStore& m)
    {
    m.makeNewData<ITDiag<Cplx>>(d.length,f.z);
    }
template void doTask(const FillCplx& f, const ITDiag<Real>& d, ManageStore& m);
template void doTask(const FillCplx& f, const ITDiag<Cplx>& d, ManageStore& m);

template<typename T>
void
doTask(const MultReal& m, ITDiag<T>& d)
    {
    d.val *= m.r;
    //use BLAS algorithm?
    for(auto& elt : d.store) elt *= m.r;
    }
template void doTask(const MultReal& m, ITDiag<Real>& d);
template void doTask(const MultReal& m, ITDiag<Cplx>& d);

template<typename T>
Real
doTask(NormNoScale, const ITDiag<T>& d)
    {
    if(d.allSame()) return std::sqrt(std::norm(d.val))*std::sqrt(d.length);
    Real nrm = 0;
    for(auto& elt : d.store) 
        nrm += std::norm(elt); //conj(elt)*elt
    return std::sqrt(nrm);
    }
template Real doTask(NormNoScale, const ITDiag<Real>& d);
template Real doTask(NormNoScale, const ITDiag<Cplx>& d);

void
doTask(Conj, ITDiag<Cplx>& d) 
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
doTask(Conj,const ITDiag<Real>& d) { }

void
doTask(TakeReal,const ITDiag<Cplx>& d, ManageStore& m) 
    { 
    if(d.allSame()) m.makeNewData<ITDiag<Real>>(d.length,d.val.real());
    else            
        {
        auto nd = m.makeNewData<ITDiag<Real>>(d.size());
        for(auto i : index(d.store)) nd->store[i] = d.store[i].real();
        }
    }

void
doTask(TakeReal, const ITDiag<Real>& ) { }

void
doTask(TakeImag,const ITDiag<Cplx>& d, ManageStore& m) 
    { 
    if(d.allSame()) m.makeNewData<ITDiag<Real>>(d.length,d.val.imag());
    else            
        {
        auto nd = m.makeNewData<ITDiag<Real>>(d.size());
        for(auto i : index(d.store)) nd->store[i] = d.store[i].imag();
        }
    }

template<typename T>
void
doTask(PrintIT<Index>& P, const ITDiag<T>& d)
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
template void doTask(PrintIT<Index>& P, const ITDiag<Real>& d);
template void doTask(PrintIT<Index>& P, const ITDiag<Cplx>& d);

bool
doTask(CheckComplex,const ITDiag<Real>& d) { return false; }

bool
doTask(CheckComplex,const ITDiag<Cplx>& d) { return true; }

template <class T>
Cplx
doTask(SumEls<Index> S, const ITDiag<T>& d) 
    { 
    if(d.allSame()) return Real(minM(S.is))*d.val;
    T sum = 0;
    for(const auto& elt : d.store)
        sum += elt;
    return sum;
    }
template Cplx doTask(SumEls<Index> S, const ITDiag<Real>& d);
template Cplx doTask(SumEls<Index> S, const ITDiag<Cplx>& d);

void
doTask(Write& W, const ITDiag<Real>& d)
    { 
    W.writeType(StorageType::ITDiagReal,d); 
    }

void
doTask(Write& W, const ITDiag<Cplx>& d)
    { 
    W.writeType(StorageType::ITDiagCplx,d); 
    }

} //namespace itensor
