//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itdata/scalar.h"
#include "itensor/util/readwrite.h"

namespace itensor {

const char*
typeNameOf(ScalarReal const& d) { return "ScalarReal"; }
const char*
typeNameOf(ScalarCplx const& d) { return "ScalarCplx"; }

} //namespace itensor
