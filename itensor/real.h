#ifndef __ITENSOR_REAL_H
#define __ITENSOR_REAL_H
#include "matrix.h"

const Real Pi = M_PI;
const Real Sqrt2 = sqrt(2);
const Real ISqrt2 = 1.0/sqrt(2);
const Real Sqrt3 = sqrt(3);
const Real ISqrt3 = 1.0/sqrt(3);

inline Real sqr(Real x) { return x*x; }

struct ApproxReal
{
    Real r;
    ApproxReal() : r(0) {}
    ApproxReal(Real _r) : r(_r) {}

    friend inline bool operator==(const ApproxReal &a,const ApproxReal &b)
    { return fabs(a.r-b.r) < 1E-12; }
    friend inline bool operator<(const ApproxReal &a,const ApproxReal &b)
    { return b.r-a.r > 1E-12; }
};



#endif
