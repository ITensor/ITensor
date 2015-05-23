//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_PRINTING_H
#define __ITENSOR_PRINTING_H

#include "itensor/util/print.h"

namespace itensor {
namespace detail {

void inline
printVal(std::ostream& s,
         double val)
    {
    if(std::fabs(val) > 1E-10)
        s << val << "\n";
    else
        s << format("%.8E\n",val);
    }

void inline
printVal(std::ostream& s,
         const std::complex<double>& val)
    {
    if(std::norm(val) > 1E-10)
        {
        auto sgn = (val.imag() < 0 ? '-' : '+');
        s << val.real() << sgn << std::fabs(val.imag()) << "i\n";
        }
    else
        {
        s << format("%.8E\n",val);
        }
    }


} //namespace detail
} //namespace itensor

#endif
