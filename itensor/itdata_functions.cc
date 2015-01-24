//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itdata_functions.h"
#include "detail/gcounter.h"
#include "detail/printing.h"
#include "lapack_wrap.h"

using std::vector;
using std::move;

namespace itensor {

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
    return move(nd);
    }

ITResult FillReal::
operator()(ITDiag<Real>& d) const
    {
    return make_newdata<ITDiag<Real>>(r_);
    }

ITResult FillReal::
operator()(const ITDiag<Complex>& d) const
    {
    return make_newdata<ITDiag<Real>>(r_);
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
    return move(nd);
    }

ITResult MultComplex::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= z_;
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
    if(a1.allSame() || a2.allSame()) Error("ITDiag plusEq allSame case not implemented");
    plusEqData(fac_,a1.data.data(),a2.data.data(),a1.data.size());
    return ITResult();
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
        auto val = scalefac*d.data[ind(is_,gc.i)];
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            detail::printVal(s_,val);
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
    auto allsame = d.allSame();
    s_ << " Diag" << (allsame ? "(all same)" : "") << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto size = minM(is_);
    for(size_t i = 0; i < size; ++i)
        {
        auto val = scalefac*(allsame ? d.val : d.data[i]);
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(int j = 1; j < is_.size(); ++j)
                {
                s_ << (1+i) << ",";
                }
            s_ << (1+i) << ") ";
            detail::printVal(s_,val);
            }
        }
    return ITResult();
    }
template ITResult PrintIT::operator()(const ITDiag<Real>& d) const;
template ITResult PrintIT::operator()(const ITDiag<Complex>& d) const;

}; //namespace itensor
