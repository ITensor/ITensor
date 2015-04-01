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

ITResult MultComplex::
operator()(const ITDense<Real>& d) const
    {
    auto nd = make_newdata<ITDense<Complex>>(d.data.begin(),d.data.end());
    operator()(*nd);
    return move(nd);
    }

void MultComplex::
operator()(ITDense<Complex>& d) const
    {
    //TODO: use BLAS algorithm
    for(auto& elt : d.data)
        elt *= z_;
    }


template<typename T>
void PrintIT::
operator()(const ITDense<T>& d) const
    {
    s_ << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = is_.r();
    if(rank == 0) 
        {
        s_ << "  ";
        detail::printVal(s_,scalefac*d.data.front());
        }

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
    }
template void PrintIT::operator()(const ITDense<Real>& d) const;
template void PrintIT::operator()(const ITDense<Complex>& d) const;

template<typename T>
void PrintIT::
operator()(const ITDiag<T>& d) const
    {
    auto allsame = d.allSame();
    s_ << " Diag" << (allsame ? "(all same)" : "") << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    if(is_.r() == 0) 
        {
        s_ << "  ";
        detail::printVal(s_,scalefac*(d.data.empty() ? d.val : d.data.front()));
        }

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
    }
template void PrintIT::operator()(const ITDiag<Real>& d) const;
template void PrintIT::operator()(const ITDiag<Complex>& d) const;

}; //namespace itensor
