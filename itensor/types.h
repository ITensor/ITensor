//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TYPES_H_
#define __ITENSOR_TYPES_H_

#include <limits>
#include <complex>

#ifndef NAN
#define NAN (std::numeric_limits<Real>::quiet_NaN())
#endif

namespace itensor {

using Real = double;
using Cplx = std::complex<Real>;
using Complex = std::complex<Real>;

}

#endif
