//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <algorithm>
#include <utility>
#include "itensor/spectrum.h"

using std::move;

namespace itensor {

struct OrderSecond
    {
    template<typename PairType>
    bool
    operator()(PairType const& i, PairType const& j) const 
        { 
        return i.second > j.second; 
        }
    };

Spectrum::
Spectrum(Args const& args) 
  : truncerr_(NAN)
    { 
    computeTruncerr(args);
    }

Spectrum::
Spectrum(Vector && eigs, Args const& args)
  : eigs_(move(eigs))
    {
    computeTruncerr(args);
    }

Spectrum::
Spectrum(Vector    && eigs, 
         QNStorage && qns,
         Args      const& args)
  : eigs_(move(eigs)),
    qns_(move(qns))
    {
    computeTruncerr(args);
    }

QN Spectrum::
qn(int n) const
    {
    if(qns_.empty()) return QN();
    return qns_.at(n-1);
    }


void Spectrum::
read(std::istream& s)
    {
    Error("Spectrum::read not currently implemented");
    //s.read((char*)&truncerr_,sizeof(truncerr_));
    //eigs_.read(s);
    //size_t sz = 0;
    //s.read((char*)&sz,sizeof(sz));
    //qns_.resize(sz);
    //for(size_t j = 0; j < sz; ++j)
    //    {
    //    s.read((char*)&qns_[j],sizeof(qns_[j]));
    //    }
    }

void Spectrum::
write(std::ostream& s) const
    {
    Error("Spectrum::read not currently implemented");
    //s.write((char*)&truncerr_,sizeof(truncerr_));
    //eigs_.write(s);
    //size_t sz = qns_.size();
    //s.write((char*)&sz,sizeof(sz));
    //for(size_t j = 0; j < sz; ++j)
    //    {
    //    s.write((char*)&qns_[j],sizeof(qns_[j]));
    //    }
    }

void Spectrum::
computeTruncerr(Args const& args)
    {
    if(args.defined("Truncerr"))
        {
        truncerr_ = args.getReal("Truncerr");
        return;
        }
    if(eigs_.size() > 0)
        {
        // Note: This only makes sense if the spectrum was normalized before truncation.
        // Otherwise this will mostly return zero.
        truncerr_ = 1.-sumels(eigs_);
        if(truncerr_ < 0) truncerr_ = 0;
        }
    }

std::ostream& 
operator<<(std::ostream & s, Spectrum const& spec)
    {
    auto& eigs = spec.eigsKept();
    auto N = eigs.size();
    if(N > 0)
        {
        decltype(N) max_show = 20;
        auto stop = std::min(N,max_show);
        s << "  Eigs kept:";
        for(auto j : range(stop))
            {
            s << format(eigs(j) > 1E-3 ? (" %.3f") : (" %.3E"), eigs(j));
            }
        s << format("\n  Trunc. error = %.3E\n", spec.truncerr());
        }
    return s;
    }


} //namespace itensor
