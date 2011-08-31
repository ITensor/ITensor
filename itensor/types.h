#ifndef __TYPES_H
#define __TYPES_H

//#define NMAX 8
static const int NMAX = 8;
static const Real MAX_CUT = 1E-15;
static const int MAX_M = 5000;

const Real Pi = M_PI;
const Real Sqrt2 = sqrt(2);
const Real ISqrt2 = 1.0/sqrt(2);
const Real Sqrt3 = sqrt(3);
const Real ISqrt3 = 1.0/sqrt(3);

inline Real sqr(Real x) { return x*x; }

//----------------------------------
//For bounds checking - can remove once implementation is tested
//#define ITENSOR_USE_AT

#ifdef  ITENSOR_USE_AT
#define GET(container,j) (container.at(j))
#else
#define GET(container,j) (container[j])
#endif
//----------------------------------

class ApproxReal
{
public:
    Real r;
    ApproxReal(Real _r = 0.0) : r(_r) {}

    friend inline bool operator==(const ApproxReal &a,const ApproxReal &b)
    { return fabs(a.r-b.r) < 1.0e-12; }
    friend inline bool operator<(const ApproxReal &a,const ApproxReal &b)
    { return b.r-a.r > 1.0e-12; }
};

#endif
