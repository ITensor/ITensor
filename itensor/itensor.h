//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "itensor/itensor_interface.h"
#include "itensor/tensor/mat.h"
#include "itensor/tensor/slicerange.h"

namespace itensor {

//
// ITensor
//
// For the ITensor class interface, see itensor_interface.h
// For the available operators see below
//

// Contract with IndexVal
// If iv = (J,n), Index J is fixed to it's nth
// value and rank decreases by 1
// (similar to summing against a Kronecker
// delta tensor \delta_{J,n})

} //namespace itensor

//See file itensor_impl.h for template/inline method implementations
#include "itensor/itensor_impl.h"


#endif
