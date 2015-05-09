//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "spectrum.h"
#include <algorithm>
#include <utility>

using std::pair;
using std::make_pair;


namespace itensor {

struct OrderSecond
    {
    using value_type = pair<QN,Real>;

    bool
    operator()(const value_type& i, const value_type& j) const 
        { 
        return i.second > j.second; 
        }
    };

Spectrum::
Spectrum(const Args& args) 
    :
    truncerr_(NAN)
    { 
    computeTruncerr(args);
    }

Spectrum::
Spectrum(const ITensor& D, const Args& args)
    :
    truncerr_(0)
    {
    Error("Spectrum ITensor constructor not yet implemented");
    //if(D.type() != ITensor::Diag) Error("Spectrum may only be constructed from Diag type ITensor.");
    //eigs_ = D.diag();
    //for(auto n = 1l; n <= eigs_.size(); ++n)
    //    eigs_(n) = sqr(eigs_(n));
    //computeTruncerr(args);
    }

Spectrum::
Spectrum(const IQTensor& D, const Args& args)
    :
    truncerr_(0)
    {
    Error("Spectrum IQTensor constructor not yet implemented");
    //std::vector<OrderSecond::value_type> eigs;
    //eigs.reserve(D.indices().front().m());

    //for(const ITensor& t : D.blocks())
    //    {
    //	if(t.type() != ITensor::Diag)
    //		Error("Spectrum may only be constructed from IQTensor containing only Diag type ITensor.");
    //    const Vec svals = t.diag();
    //    const QN q = itensor::qn(D,t.indices().front());
    //    for(int n = 1; n <= svals.size(); ++n)
    //        {
    //        eigs.push_back(std::make_pair(q,sqr(svals(n))));
    //        }
    //    }
    //std::sort(eigs.begin(),eigs.end(),OrderSecond());

    //qns_.resize(eigs.size());
    //eigs_.ReDimension(eigs.size());
    //for(size_t j = 0; j < eigs.size(); ++j)
    //    {
    //    qns_.at(j) = eigs.at(j).first;
    //    eigs_[j] = eigs.at(j).second;
    //    }
    //computeTruncerr(args);
    }

Spectrum::
Spectrum(const Vec& eigs, const Args& args)
    :
    eigs_(eigs)
    {
    computeTruncerr(args);
    }


Spectrum::
Spectrum(const Vec& eigs, 
         const QNStorage& qns,
         const Args& args)
    :
    eigs_(eigs),
    qns_(qns)
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
computeTruncerr(const Args& args)
    {
    if(args.defined("Truncerr"))
        {
        truncerr_ = args.getReal("Truncerr");
        }
    else
        {
        if(eigs_.size() > 0)
            {
        	// Note: This only makes sense if the spectrum was normalized before truncation.
        	// Otherwise this will mostly return zero.
            truncerr_ = 1.-sumels(eigs_);
            if(truncerr_ < 0) truncerr_ = 0;
            }
        }
    }

std::ostream& 
operator<<(std::ostream & s, const Spectrum& spec)
    {
    const Vec& eigs = spec.eigsKept();
    auto N = eigs.size();
    if(N > 0)
        {
        long max_show = 20;
        auto stop = std::min(N,max_show);
        s << "  Eigs kept: ";
        for(auto j = 1l; j <= stop; ++j)
            {
            s << format(eigs(j) > 1E-3 ? ("%.3f") : ("%.3E"), eigs(j));
            s << ((j != stop) ? ", " : "\n");
            }
        s << format("  Trunc. error = %.3E\n", spec.truncerr());
        }
    return s;
    }



} //namespace itensor
