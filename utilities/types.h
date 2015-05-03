//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_TYPES_H_
#define __ITENSOR_TYPES_H_

#include <complex>
#include <limits>

#define Foreach(X,Y) for(X : Y)

namespace itensor {

using Real = double;
using Complex = std::complex<Real>;
using Cplx = std::complex<Real>;

#ifndef NAN
#define NAN (std::numeric_limits<Real>::quiet_NaN())
#endif

};

#include <array>
namespace itensor {
using std::array;
};

#include <memory>
namespace itensor {
using std::shared_ptr;
using std::make_shared;
};

#include <random>
namespace itensor {
using std::mt19937;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
};

#include <functional>
namespace itensor {
using std::function;
};

#endif
