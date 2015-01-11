//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itdata_functions.h"
#include "detail/gcounter.h"
#include "lapack_wrap.h"

using std::vector;

namespace itensor {

ITResult Contract::
operator()(const ITDense<Real>& a1,
           const ITDense<Real>& a2)
    {
    const auto& Lis = *Lis_;
    const auto& Ris = *Ris_;
    const auto& Lind = *Lind_;
    const auto& Rind = *Rind_;

    long ncont = 0;
    for(const auto& i : Lind) if(i < 0) ++ncont;
    long nuniq = Lis.r()+Ris.r()-2*ncont;
    vector<Index> newind(nuniq);

    long nn = 0;
    for(int j = 0; j < Lis.r(); ++j)
        {
        if(Lind[j] > 0) 
            {
            newind[nn++] = Lis[j];
            }
        }
    for(int j = 0; j < Ris.r(); ++j)
        {
        if(Rind[j] > 0) 
            {
            newind[nn++] = Ris[j];
            }
        }
    auto comp = [](const Index& i1, const Index& i2) { return i1 > i2; };
    std::sort(newind.begin(),newind.end(),comp);
    Nis_ = IndexSet(std::move(newind));
    
    Label Nind(nuniq);
    for(size_t i = 0; i < Nis_.r(); ++i)
        {
        auto j = findindex(*Lis_,Nis_[i]);
        if(j >= 0)
            {
            Nind[i] = (*Lind_)[j];
            }
        else
            {
            j = findindex(*Ris_,Nis_[i]);
            Nind[i] = (*Rind_)[j];
            }
        }

    //PRI(Lind);
    //PRI(Rind);
    //PRI(Nind);

    auto res = make_newdata<ITDense<Real>>(area(Nis_),0.);
    auto t1 = make_tensorref(a1.data.data(),Lis),
         t2 = make_tensorref(a2.data.data(),Ris),
         tr = make_tensorref(res->data.data(),Nis_);
    contractloop(t1,Lind,t2,Rind,tr,Nind);
    return std::move(res);
    }

ITResult FillReal::
operator()(ITDense<Real>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),r_);
    return ITResult();
    }

ITResult FillReal::
operator()(const ITDense<Complex>& d) const
    {
    auto nd = make_newdata<ITDense<Real>>(d.data.size());
    operator()(*nd);
    return std::move(nd);
    }

ITResult FillReal::
operator()(ITDiag<Real>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),r_);
    return ITResult();
    }

ITResult FillReal::
operator()(const ITDiag<Complex>& d) const
    {
    return make_newdata<ITDiag<Real>>(d.data.size(),r_);
    }

ITResult FillCplx::
operator()(const ITDense<Real>& d) const
    {
    return make_newdata<ITDense<Complex>>(d.data.size(),z_);
    }

ITResult FillCplx::
operator()(ITDense<Complex>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),z_);
    return ITResult();
    }

ITResult MultComplex::
operator()(const ITDense<Real>& d) const
    {
    auto nd = make_newdata<ITDense<Complex>>(d.data.begin(),d.data.end());
    operator()(*nd);
    return std::move(nd);
    }

ITResult MultComplex::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= z_;
    return ITResult();
    }

ITResult MultReal::
operator()(ITDense<Real>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= r_;
    return ITResult();
    }

ITResult MultReal::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= r_;
    return ITResult();
    }

void
plusEqData(Real fac, Real *d1, const Real *d2, LAPACK_INT size)
    {
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,d2,&inc,d1,&inc);
    }

ITResult PlusEQ::
operator()(ITDense<Real>& a1,
           const ITDense<Real>& a2)
    {
#ifdef DEBUG
    if(a1.data.size() != a2.data.size()) Error("Mismatched sizes in plusEq");
#endif
    if(permute_)
        {
        auto ref1 = tensorref<Real,IndexSet>(a1.data.data(),*is1_),
             ref2 = tensorref<Real,IndexSet>(a2.data.data(),*is2_);
        auto f = fac_;
        auto add = [f](Real& r1, Real r2) { r1 += f*r2; };
        reshape(ref2,*P_,ref1,add);
        }
    else
        {
        plusEqData(fac_,a1.data.data(),a2.data.data(),a1.data.size());
        }
    return ITResult();
    }

ITResult PlusEQ::
operator()(ITDiag<Real>& a1,
           const ITDiag<Real>& a2)
    {
#ifdef DEBUG
    if(a1.data.size() != a2.data.size()) Error("Mismatched sizes in plusEq");
#endif
    plusEqData(fac_,a1.data.data(),a2.data.data(),a1.data.size());
    return ITResult();
    }

void
printVal(std::ostream& s,
         Real val)
    {
    if(std::fabs(val) > 1E-10)
        s << val << "\n";
    else
        s << format("%.8E\n",val);
    }

void
printVal(std::ostream& s,
         const Complex& val)
    {
    if(std::norm(val) > 1E-10)
        {
        auto sgn = (val.imag() < 0 ? '-' : '+');
        s << val.real() << sgn << std::fabs(val.imag()) << "i\n";
        }
    else
        {
        s << format("%.8E\n",val);
        }
    }

template<typename T>
ITResult PrintIT::
operator()(const ITDense<T>& d) const
    {
    s_ << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) return ITResult();

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,is_.dim(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = d.data[ind(is_,gc.i)];
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            printVal(s_,val);
            }
        }
    return ITResult();
    }
template ITResult PrintIT::operator()(const ITDense<Real>& d) const;
template ITResult PrintIT::operator()(const ITDense<Complex>& d) const;

template<typename T>
ITResult PrintIT::
operator()(const ITDiag<T>& d) const
    {
    s_ << " Diag}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    for(size_t i = 0; i < d.data.size(); ++i)
        {
        auto val = d.data[i];
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(int j = 1; j < is_.size(); ++j)
                {
                s_ << (1+i) << ",";
                }
            s_ << (1+i) << ") ";
            printVal(s_,val);
            }
        }
    return ITResult();
    }
template ITResult PrintIT::operator()(const ITDiag<Real>& d) const;
template ITResult PrintIT::operator()(const ITDiag<Complex>& d) const;

}; //namespace itensor
