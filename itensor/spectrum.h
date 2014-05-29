//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPECTRUM_H
#define __ITENSOR_SPECTRUM_H

#include "iqtensor.h"

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

    typedef std::vector<QN>
    QNStorage;

    Spectrum(const OptSet& opts = Global::opts());

    Spectrum(const ITensor& D, const OptSet& opts = Global::opts());

    Spectrum(const IQTensor& D, const OptSet& opts = Global::opts());

    Spectrum(const Vector& eigs, const OptSet& opts = Global::opts());

    Spectrum(const Vector& eigs, 
             const QNStorage& qns,
             const OptSet& opts = Global::opts());

    QN
    qn(int n) const;

    Real
    eig(int n) const { return eigsKept_(n); }

    Real 
    truncerr() const { return truncerr_; }

    int
    numEigsKept() const { return eigsKept_.Length(); }

    bool
    hasQNs() const { return !qns_.empty(); }

    const Vector& 
    eigsKept() const { return eigsKept_; }

    const QNStorage&
    qn() const { return qns_; }


    //
    // Other Methods
    //

    void 
    truncerr(Real val) { truncerr_ = val; }

    void 
    eigsKept(const Vector& val) { eigsKept_ = val; }

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:

    /////////////////

    Vector eigsKept_;
    Real truncerr_;
    std::vector<QN> qns_;

    /////////////////

    void
    computeTruncerr(const OptSet& opts);

    }; //class Spectrum

std::ostream& 
operator<<(std::ostream & s, const Spectrum& spec);

}; //namespace itensor


#endif
