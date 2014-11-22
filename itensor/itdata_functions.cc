#include "itdata_functions.h"
#include "detail/gcounter.h"

namespace itensor {

NewData FillReal::
operator()(ITDense<Real>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),r_);
    return NewData();
    }

NewData FillReal::
operator()(const ITDense<Complex>& d) const
    {
    auto nd = make_newdata<ITDense<Real>>(d.data.inds());
    operator()(*nd);
    return std::move(nd);
    }

NewData FillCplx::
operator()(const ITDense<Real>& d) const
    {
    auto nd = make_newdata<ITDense<Complex>>(d.data.inds());
    operator()(*nd);
    return std::move(nd);
    }

NewData FillCplx::
operator()(ITDense<Complex>& d) const
    {
    std::fill(d.data.begin(),d.data.end(),z_);
    return NewData();
    }

//NewData MultComplex::
//operator()(const ITDense<Real>& d) const
//    {
//    auto nd = new ITDense<Complex>(d);
//    btas::scal(z_,nd->t_);
//    return NewData(nd);
//    }
//
//NewData MultComplex::
//operator()(ITDense<Complex>& d) const
//    {
//    btas::scal(z_,d.t_);
//    return NewData();
//    }
//
//NewData MultComplex::
//operator()(const ITScalar<Real>& d) const
//    {
//    auto nd = new ITScalar<Complex>(z_*d.x_);
//    return NewData(nd);
//    }
//
//NewData MultComplex::
//operator()(ITScalar<Complex>& d) const
//    {
//    d.x_ *= z_;
//    return NewData();
//    }
//
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

NewData PrintIT::
operator()(const ITDense<Real>& d) const
    {
    s_ << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = d.data.r();
    if(rank == 0) return NewData();

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,d.data.n(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = d.data(gc.i);
        if(fabs(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            if(fabs(val) > 1E-10)
                s_ << val << "\n";
            else
                s_ << format("%.8E\n",val);
            }
        }
    return NewData();
    }

NewData PrintIT::
operator()(const ITDense<Complex>& d) const 
    {
    s_ << "}\n";
    Real scalefac = 1.0;
    if(!x_.isTooBigForReal()) scalefac = x_.real0();
    else s_ << "  (omitting too large scale factor)\n";

    auto rank = d.data.r();
    if(rank == 0) return NewData();

    auto gc = detail::GCounter(0,rank-1,0);
    for(int i = 0; i < rank; ++i)
        gc.setInd(i,0,d.data.n(i)-1);

    for(; gc.notDone(); ++gc)
        {
        auto val = d.data(gc.i)*scalefac;
        if(std::norm(val) > Global::printScale())
            {
            s_ << "  (";
            for(auto ii = gc.i.mini(); ii <= gc.i.maxi(); ++ii)
                {
                s_ << (1+gc.i(ii));
                if(ii < gc.i.maxi()) s_ << ",";
                }
            s_ << ") ";

            if(std::norm(val) > 1E-10)
                {
                auto sgn = (val.imag() < 0 ? '-' : '+');
                s_ << val.real() << sgn << fabs(val.imag()) << "i\n";
                }
            else
                {
                s_ << format("%.8E\n",val);
                }
            }
        }
    return NewData();
    }



}; //namespace itensor
