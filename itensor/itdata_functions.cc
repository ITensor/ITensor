//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itdata_functions.h"
#include "detail/gcounter.h"
#include "lapack_wrap.h"
#include "contract.h"

namespace itensor {

NewData Contract::
operator()(const ITDense<Real>& a1,
           const ITDense<Real>& a2) const
    {
    auto res = make_newdata<ITDense<Real>>(area(Nis_),0.);
    auto t1 = make_tensorref(a1.data.data(),Lis_),
         t2 = make_tensorref(a2.data.data(),Ris_),
         tr = make_tensorref(res->data.data(),Nis_);
    contractloop(t1,Lind_,t2,Rind_,tr,Nind_);
    return std::move(res);
    }

NewData Contract::
operator()(const ITDense<Real>& d,
           const ITCombiner& C) const
    {
    Error("Not implemented");
    return NewData();
    }

NewData FillReal::
operator()(ITDense<Real>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),r_);
    return NewData();
    }

NewData FillReal::
operator()(const ITDense<Complex>& d) const
    {
    auto nd = make_newdata<ITDense<Real>>(d.data.size());
    operator()(*nd);
    return std::move(nd);
    }

NewData FillReal::
operator()(ITDiag<Real>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),r_);
    return NewData();
    }

NewData FillReal::
operator()(const ITDiag<Complex>& d) const
    {
    return make_newdata<ITDiag<Real>>(d.data.size(),r_);
    }

NewData FillCplx::
operator()(const ITDense<Real>& d) const
    {
    return make_newdata<ITDense<Complex>>(d.data.size(),z_);
    }

NewData FillCplx::
operator()(ITDense<Complex>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),z_);
    return NewData();
    }

NewData MultComplex::
operator()(const ITDense<Real>& d) const
    {
    auto nd = make_newdata<ITDense<Complex>>(d.data.begin(),d.data.end());
    operator()(*nd);
    return std::move(nd);
    }

NewData MultComplex::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= z_;
    return NewData();
    }

NewData MultReal::
operator()(ITDense<Real>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= r_;
    return NewData();
    }

NewData MultReal::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= r_;
    return NewData();
    }

void
plusEqData(Real fac, Real *d1, const Real *d2, LAPACK_INT size)
    {
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,d2,&inc,d1,&inc);
    }

NewData PlusEQ::
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
    return NewData();
    }

NewData PlusEQ::
operator()(ITDiag<Real>& a1,
           const ITDiag<Real>& a2)
    {
#ifdef DEBUG
    if(a1.data.size() != a2.data.size()) Error("Mismatched sizes in plusEq");
#endif
    plusEqData(fac_,a1.data.data(),a2.data.data(),a1.data.size());
    return NewData();
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
NewData PrintIT::
operator()(const ITDense<T>& d) const
    {
    s_ << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) return NewData();

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
    return NewData();
    }
template NewData PrintIT::operator()(const ITDense<Real>& d) const;
template NewData PrintIT::operator()(const ITDense<Complex>& d) const;

template<typename T>
NewData PrintIT::
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
    return NewData();
    }
template NewData PrintIT::operator()(const ITDiag<Real>& d) const;
template NewData PrintIT::operator()(const ITDiag<Complex>& d) const;

}; //namespace itensor
