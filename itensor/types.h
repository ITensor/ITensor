//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TYPES_H_
#define __ITENSOR_TYPES_H_

#include <limits>
#include "util/cplx_literal.h"

#ifndef NAN
#define NAN (std::numeric_limits<Real>::quiet_NaN())
#endif

namespace itensor {

using Real = double;
using Cplx = std::complex<double>;
using Complex = std::complex<double>;

static const Cplx Complex_1 = Cplx(1,0);
static const Cplx Complex_i = Cplx(0,1);
static const Cplx Cplx_1 = Cplx(1,0);
static const Cplx Cplx_i = Cplx(0,1);

}

#endif
