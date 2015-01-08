//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "simpletensor.h"
#include "lapack_wrap.h"

namespace itensor {

void
plusEqData(double fac,
           const double *d1,
           double *d2,
           LAPACK_INT size)
    {
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,d1,&inc,d2,&inc);
    }

};
