//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPECTRUM_H
#define __ITENSOR_SPECTRUM_H

#include "itensor/tensor/vec.h"
#include "itensor/iqtensor.h"

namespace itensor {

//
// Spectrum
//
// Stores density matrix eigenvalue spectrum following a call to
// svd, denmatDecomp, etc.
//

class Spectrum
    {
    public:
    using QNStorage = std::vector<QN>;
    private:
    Vector eigs_;
    Real truncerr_;
    QNStorage qns_;
    public:

    Spectrum(Args const& args = Args::global());

    Spectrum(Vector && eigs, Args const& args = Args::global());

    Spectrum(Vector && eigs, 
             QNStorage && qns,
             Args const& args = Args::global());

    Real 
    truncerr() const { return truncerr_; }

    Vector const&
    eigs() const { return eigs_; }

    Vector const& 
    eigsKept() const { return eigs(); }

    //1-indexed
    Real
    eig(int n) const { return eigs_(n-1); }

    int
    numEigsKept() const { return eigs_.size(); }

    bool
    hasQNs() const { return !qns_.empty(); }

    //1-indexed
    QN
    qn(int n) const;

    QNStorage const&
    qns() const { return qns_; }

    int
    size() const { return eigs_.size(); }

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:
     
    void
    computeTruncerr(Args const& args);

    }; //class Spectrum

std::ostream& 
operator<<(std::ostream & s,Spectrum const& spec);

} //namespace itensor


#endif
