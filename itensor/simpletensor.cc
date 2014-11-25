//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "simpletensor.h"
#include "lapack_wrap.h"

namespace itensor {

void
plusEq(Real fac,
       const simpletensor<Real>& t1,
       simpletensor<Real>& t2)
    {
#ifdef DEBUG
    if(t1.size() != t2.size()) Error("Mismatched sizes in plusEq");
#endif
    LAPACK_INT size = t1.size();
    LAPACK_INT inc = 1;
    daxpy_wrapper(&size,&fac,t1.data(),&inc,t2.data(),&inc);
    }

};
