//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "spectrum.h"
#include <algorithm>
#include <utility>

struct OrderSecond
    {
    using value_type = pair<QN,Real>;
    bool
    operator()(const value_type& i, const value_type& j) const 
        { 
        return i.second > j.second; 
        }
    };

namespace itensor {

Spectrum::
Spectrum(const OptSet& opts) 
    :
    truncerr_(NAN)
    { 
    }

Spectrum::
Spectrum(const ITensor& D, const OptSet& opts)
    :
    truncerr_(0)
    {
    if(D.type() != ITensor::Diag)
        Error("Spectrum may only be constructed from Diag type ITensor.");

    eigsKept_ = D.diag();
    truncerr_ = 1.-eigsKept_.sumels();
    if(truncerr < 0) truncerr_ = 0;
    }

Spectrum::
Spectrum(const IQTensor& D, const OptSet& opts)
    :
    truncerr_(0)
    {
    std::vector<OrderSecond::value_type> eigs;
    eigs.reserve(D.indices().front().m());

    Foreach(const ITensor& t, D.blocks())
        {
        const Vector teigs = t.diag();
        const QN q = qn(D,t.indices().front());
        for(int n = 1; n <= teigs.Length(); ++n)
            {
            eigs.push_back(std::make_pair(q,teigs(n)));
            }
        }

    std::sort(eigs.begin(),eigs.end(),OrderSecond());
    }

Spectrum::
Spectrum(const Vector& eigs, const OptSet& opts)
    {
    }

Spectrum::
Spectrum(const Vector& eigs, 
         const QNStorage& qns,
         const OptSet& opts)
    {
    }

void Spectrum::
read(std::istream& s)
    {
    s.read((char*)&truncerr_,sizeof(truncerr_));
    eigsKept_.read(s);
    }

void Spectrum::
write(std::ostream& s) const
    {
    s.write((char*)&truncerr_,sizeof(truncerr_));
    eigsKept_.write(s);
    }

std::ostream& 
operator<<(std::ostream & s, const Spectrum& spec)
    {
    const Vector& eigs = spec.eigsKept();
    const int N = eigs.Length();
    if(N > 0)
        {
        const int max_show = 20;
        int stop = min(N,max_show);
        s << "  Eigs kept: ";
        for(int j = 1; j <= stop; ++j)
            {
            s << format(eigs(j) > 1E-3 ? ("%.3f") : ("%.3E"), eigs(j));
            s << ((j != stop) ? ", " : "\n");
            }
        s << format("  Trunc. error = %.3E\n", spec.truncerr());
        }
    return s;
    }



}; //namespace itensor
