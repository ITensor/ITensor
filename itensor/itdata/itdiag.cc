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

void
diagDense(const ITDiag<Real>& d,
          const IndexSet& dis,
          const Label& dind,
          const ITReal& t,
          const IndexSet& tis,
          const Label& tind,
          IndexSet& Nis,
          ManageStore& m,
          bool RealonLeft)
    {
    //if(area(tis)==1) //Dense tensor is a scalar
    //    {
    //    scalefac_ = R[0];
    //    if(RealOnLeft) assignPointerRtoL();
    //    return;
    //    }
    //if(area(dis)==1) //Diag tensor is a scalar
    //    {
    //    Error("Not implemented");
    //    return;
    //    }

    long t_cstride = 0; //total t-stride of contracted inds of t
    size_t ntu = 0; //number uncontracted inds of t
    assert(tind.size() == tis.size());
    for(auto j : index(tind))
        {
        //if index j is contracted, add its stride to t_cstride:
        if(tind[j] < 0) t_cstride += tis.stride(j);
        else            ++ntu;
        }

    long d_ustride = 0; //total result-stride of uncontracted inds of d
    for(auto i : index(Nis))
        {
        auto j = findindex(dis,Nis[i]);
        if(j >= 0) d_ustride += Nis.stride(i);
        }


    if(ntu > 0)
        {
        vector<long> tstride(ntu,0),
                     rstride(ntu,0);
        detail::GCounter C(0,ntu,0);
        size_t n = 0;
        for(auto j : index(tind))
            {
            if(tind[j] > 0)
                {
#ifdef DEBUG
                if(n >= ntu) Error("n out of range");
#endif
                C.setInd(n,0,tis.extent(j)-1);
                tstride.at(n) = tis.stride(j);
                auto k = findindex(Nis,tis[j]);
#ifdef DEBUG
                if(k < 0) Error("Index not found");
#endif
                rstride.at(n) = Nis.stride(k);
                ++n;
                }
            }
        auto nd = m.makeNewData<ITReal>(area(Nis),0.);
        auto pr = MAKE_SAFE_PTR(nd->data(),nd->size());
        auto pt = MAKE_SAFE_PTR(t.data(),t.size());

        if(d.allSame())
            {
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(auto i : count(ntu))
                    {
                    auto ii = C.i[i];
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(auto J : count(d.length))
                    {
                    pr[J*d_ustride+roffset] += d.val*pt[J*t_cstride+toffset];
                    }
                }
            }
        else
            {
            auto pd = MAKE_SAFE_PTR(d.data(),d.size());
            assert(d.size() == d.length);
            for(;C.notDone();++C)
                {
                size_t roffset = 0,
                       toffset = 0;
                for(auto i : count(ntu))
                    {
                    auto ii = C.i.fast(i);
                    toffset += ii*tstride[i];
                    roffset += ii*rstride[i];
                    }
                for(auto J : count(d.length))
                    {
                    pr[J*d_ustride+roffset] += pd[J]*pt[J*t_cstride+toffset];
                    }
                }
            }
        }
    else
        {
        //all of t's indices contracted with d
        //result will be diagonal
        if(d_ustride == 0) //all of d's inds contracted
            {
            // o scalar if all of d's inds contracted also
            Real val = 0;
            auto pt = MAKE_SAFE_PTR(t.data(),t.size());
            if(d.allSame())
                {
                for(auto J : count(d.length))
                    val += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(d.length == d.size());
                auto pd = MAKE_SAFE_PTR(d.data(),d.size());
                for(auto J : count(d.length))
                    val += pd[J]*pt[J*t_cstride];
                }
            m.makeNewData<ITDiag<Real>>(1,val);
            }
        else //some of d's inds uncontracted
            {
            // o element-wise product of d's data and t's diagonal
            auto nd = m.makeNewData<ITDiag<Real>>(d.length);
            auto pr = MAKE_SAFE_PTR(nd->data(),nd->size());
            auto pt = MAKE_SAFE_PTR(t.data(),t.size());
            if(d.allSame())
                {
                for(auto J : count(d.length))
                    pr[J] += d.val*pt[J*t_cstride];
                }
            else
                {
                assert(d.length == d.size());
                auto pd = MAKE_SAFE_PTR(d.data(),d.size());
                for(auto J : count(d.length))
                    pr[J] += pd[J]*pt[J*t_cstride];
                }
            }
        }
    }

void
doTask(Contract<Index>& C,
       const ITReal& t,
       const ITDiag<Real>& d,
       ManageStore& m)
    { 
    contractIS(C.Lis,C.Lind,C.Ris,C.Rind,C.Nis,true);
    diagDense(d,C.Ris,C.Rind,t,C.Lis,C.Lind,C.Nis,m,true);
    }
void
doTask(Contract<Index>& C,
       const ITDiag<Real>& d,
       const ITReal& t,
       ManageStore& m)
    {
    contractIS(C.Lis,C.Lind,C.Ris,C.Rind,C.Nis,true);
    diagDense(d,C.Lis,C.Lind,t,C.Ris,C.Rind,C.Nis,m,false);
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
