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
#ifndef __ITENSOR_TYPES_H_
#define __ITENSOR_TYPES_H_

#include <limits>
#include "itensor/util/stdx.h"
#include "itensor/util/cplx_literal.h"
#include "itensor/util/print.h"

#ifndef NAN
#define NAN (std::numeric_limits<Real>::quiet_NaN())
#endif

namespace itensor {

using Real = double;
using Cplx = std::complex<double>;
using Complex = std::complex<double>;

const Cplx Complex_1 = Cplx(1,0);
const Cplx Complex_i = Cplx(0,1);
const Cplx Cplx_1 = Cplx(1,0);
const Cplx Cplx_i = Cplx(0,1);

inline Real& 
realRef(Cplx & z) { return reinterpret_cast<Real*>(&z)[0]; }

inline Real& 
imagRef(Cplx & z) { return reinterpret_cast<Real*>(&z)[1]; }

void inline
applyConj(Real & r) { }

void inline
applyConj(Cplx & z) { imagRef(z) *= -1; }

std::string inline
formatVal(double val)
    {
    if(val == 0.0 || val >= 1E-7)
        {
        return format("%.7f",val);
        }
    else if(val <= -1E-7)
        {
        return format("%.6f",val);
        }
    return format("%.8E",val);
    }

std::string inline
formatVal(Cplx const& val)
    {
    auto sgn = (val.imag() < 0 ? '-' : '+');
    auto nrm = std::norm(val);
    if(nrm == 0. || nrm > 1E-10)
        {
        return format("%f%s%fi",val.real(),sgn,std::fabs(val.imag()));
        }
    else
        {
        return format("%.8E%s%.8Ei",val.real(),sgn,std::fabs(val.imag()));
        }
    }

template<typename T, class=stdx::require<std::is_same<T,Real>>>
constexpr const char* 
typeName(int=0) { return "Real"; }

template<typename T, class=stdx::require<std::is_same<T,Cplx>>>
constexpr const char* 
typeName(long=0) { return "Cplx"; }

}

#endif
