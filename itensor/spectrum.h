//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_SPECTRUM_H
#define __ITENSOR_SPECTRUM_H

#include "itensor/tensor/vec.h"
#include "itensor/qn.h"

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
