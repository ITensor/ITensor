//
// Distributed under the ITensor Library License, Version 1.1.
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
    typedef pair<QN,Real> value_type;

    bool
    operator()(const value_type& i, const value_type& j) const 
        { 
        return i.second > j.second; 
        }
    };

Spectrum::
Spectrum(const OptSet& opts) 
    :
    truncerr_(NAN)
    { 
    computeTruncerr(opts);
    }

Spectrum::
Spectrum(const ITensor& D, const OptSet& opts)
    :
    truncerr_(0)
    {
    if(D.type() != ITensor::Diag) Error("Spectrum may only be constructed from Diag type ITensor.");
    eigs_ = D.diag();
    for(int n = 1; n <= eigs_.Length(); ++n)
        eigs_(n) = sqr(eigs_(n));
    computeTruncerr(opts);
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
    	if(t.type() != ITensor::Diag)
    		Error("Spectrum may only be constructed from IQTensor containing only Diag type ITensor.");
        const Vector svals = t.diag();
        const QN q = itensor::qn(D,t.indices().front());
        for(int n = 1; n <= svals.Length(); ++n)
            {
            eigs.push_back(std::make_pair(q,sqr(svals(n))));
            }
        }
    std::sort(eigs.begin(),eigs.end(),OrderSecond());

    qns_.resize(eigs.size());
    eigs_.ReDimension(eigs.size());
    for(size_t j = 0; j < eigs.size(); ++j)
        {
        qns_.at(j) = eigs.at(j).first;
        eigs_[j] = eigs.at(j).second;
        }
    computeTruncerr(opts);
    }

Spectrum::
Spectrum(const Vector& eigs, const OptSet& opts)
    :
    eigs_(eigs)
    {
    computeTruncerr(opts);
    }


Spectrum::
Spectrum(const Vector& eigs, 
         const QNStorage& qns,
         const OptSet& opts)
    :
    eigs_(eigs),
    qns_(qns)
    {
    computeTruncerr(opts);
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
    s.read((char*)&truncerr_,sizeof(truncerr_));
    eigs_.read(s);
    size_t sz = 0;
    s.read((char*)&sz,sizeof(sz));
    qns_.resize(sz);
    for(size_t j = 0; j < sz; ++j)
        {
        s.read((char*)&qns_[j],sizeof(qns_[j]));
        }
    }

void Spectrum::
write(std::ostream& s) const
    {
    s.write((char*)&truncerr_,sizeof(truncerr_));
    eigs_.write(s);
    size_t sz = qns_.size();
    s.write((char*)&sz,sizeof(sz));
    for(size_t j = 0; j < sz; ++j)
        {
        s.write((char*)&qns_[j],sizeof(qns_[j]));
        }
    }

void Spectrum::
computeTruncerr(const OptSet& opts)
    {
    if(opts.defined("Truncerr"))
        {
        truncerr_ = opts.getReal("Truncerr");
        }
    else
        {
        if(eigs_.Length() > 0)
            {
        	// Note: This only makes sense if the spectrum was normalized before truncation.
        	// Otherwise this will mostly return zero.
            truncerr_ = 1.-eigs_.sumels();
            if(truncerr_ < 0) truncerr_ = 0;
            }
        }
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
