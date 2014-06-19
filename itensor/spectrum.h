//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SPECTRUM_H
#define __ITENSOR_SPECTRUM_H

#include "iqtensor.h"
#include "iterpair.h"

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

    const QNStorage&
    qns() const { return qns_; }

    Real
    eig(int n) const { return eigs_(n); }

    const Vector&
    eigs() const { return eigs_; }

    Real 
    truncerr() const { return truncerr_; }

    int
    size() const { return eigs_.Length(); }

    bool
    hasQNs() const { return !qns_.empty(); }


    //
    // Other Methods
    //

    const Vector& 
    eigsKept() const { return eigs_; }

    int
    numEigsKept() const { return eigs_.Length(); }

    void 
    truncerr(Real val) { truncerr_ = val; }

    void 
    eigsKept(const Vector& val) { eigs_ = val; }

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:

    /////////////////

    Vector eigs_;
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
