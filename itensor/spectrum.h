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

    Spectrum(Vector const& eigs, Args const& args = Args::global());

    Spectrum(Vector const& eigs, 
             QNStorage const& qns,
             Args const& args = Args::global());

    //Spectrum(ITensor const& D, Args const& args = Args::global());

    //Spectrum(IQTensor const& D, Args const& args = Args::global());


    //1-indexed
    QN
    qn(int n) const;

    QNStorage const&
    qns() const { return qns_; }

    //1-indexed
    Real
    eig(int n) const { return eigs_(n-1); }

    Vector const&
    eigs() const { return eigs_; }

    Real 
    truncerr() const { return truncerr_; }

    int
    size() const { return eigs_.size(); }

    bool
    hasQNs() const { return !qns_.empty(); }


    //
    // Other Methods
    //

    Vector const& 
    eigsKept() const { return eigs_; }

    int
    numEigsKept() const { return eigs_.size(); }

    void 
    truncerr(Real val) { truncerr_ = val; }

    void 
    eigsKept(Vector const& val) { eigs_ = val; }

    void 
    eigsKept(Vector&& val) { eigs_ = std::move(val); }

    void 
    read(std::istream& s);
    void 
    write(std::ostream& s) const;

    private:
     
    void
    computeTruncerr(const Args& args);

    }; //class Spectrum

std::ostream& 
operator<<(std::ostream & s,Spectrum const& spec);

} //namespace itensor


#endif
