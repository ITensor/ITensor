//$$ precisio.h                          floating point constants

#ifndef PRECISION_LIB
#define PRECISION_LIB 0

//#define USING_DOUBLE		   /* Kludge RMN */
#define USING_FLOAT		   /* Kludge SRW */
//#define SystemV
//#include <values.h>
//#include <limits>

#ifndef SystemV                    // if there is float.h
#include <float.h>


#ifdef USING_FLOAT


class FloatingPointPrecision
{
public:
   static int Dig()
      { return FLT_DIG; }        // number of decimal digits or precision
   static Real Epsilon()
      { return FLT_EPSILON; }    // smallest number such that 1+Eps!=Eps
   static int Mantissa()
      { return FLT_MANT_DIG; }   // bits in mantisa
   static Real Maximum()
      { return FLT_MAX; }        // maximum value
   static int MaximumDecimalExponent()
      { return FLT_MAX_10_EXP; } // maximum decimal exponent
   static int MaximumExponent()
      { return FLT_MAX_EXP; }    // maximum binary exponent
   static Real Minimum()
      { return FLT_MIN; }        // minimum positive value
   static int MinimumDecimalExponent()
      { return FLT_MIN_10_EXP; } // minimum decimal exponent
   static int MinimumExponent()
      { return FLT_MIN_EXP; }    // minimum binary exponent
   static int Radix()
      { return FLT_RADIX; }      // exponent radix
   static int Rounds()
      { return FLT_ROUNDS; }     // addition rounding (1 = does round)
//   FREE_CHECK(FloatingPointPrecision)
};

#endif


#ifdef USING_DOUBLE

class FloatingPointPrecision
{
public:
   static int Dig()
      { return DBL_DIG; }        // number of decimal digits or precision
   static Real Epsilon()
      { return DBL_EPSILON; }    // smallest number such that 1+Eps!=Eps
   static int Mantissa()
      { return DBL_MANT_DIG; }   // bits in mantisa
   static Real Maximum()
      { return DBL_MAX; }        // maximum value
   static int MaximumDecimalExponent()
      { return DBL_MAX_10_EXP; } // maximum decimal exponent
   static int MaximumExponent()
      { return DBL_MAX_EXP; }    // maximum binary exponent
   static Real Minimum()
   {
// #ifdef __BCPLUSPLUS__
        return 2.225074e-308;     // minimum positive value
// #else
//        return DBL_MIN;
// #endif
   }
   static int MinimumDecimalExponent()
      { return DBL_MIN_10_EXP; } // minimum decimal exponent
   static int MinimumExponent()
      { return DBL_MIN_EXP; }    // minimum binary exponent
   static int Radix()
      { return FLT_RADIX; }      // exponent radix
   static int Rounds()
      { return FLT_ROUNDS; }     // addition rounding (1 = does round)
//   FREE_CHECK(FloatingPointPrecision)
};

#endif

#endif

#ifdef SystemV                    // if there is no float.h

#ifdef USING_FLOAT

class FloatingPointPrecision
{
public:
   static Real Epsilon()
      { return pow(2.0,1-FSIGNIF); }  // smallest number such that 1+Eps!=Eps
   static Real Maximum()
      { return MAXFLOAT; }        // maximum value
   static Real Minimum()
      { return MINFLOAT; }        // minimum positive value
//   FREE_CHECK(FloatingPointPrecision)
};

#endif


#ifdef USING_DOUBLE

class FloatingPointPrecision
{
public:
   static Real Epsilon()
     // { return pow(2.0,1-DSIGNIF); }  // smallest number such that 1+Eps!=Eps
      { return 2.2204460492503131e-16; }  
			// smallest number such that 1+Eps!=Eps
   static Real Maximum()
      { return MAXDOUBLE; }          // maximum value
   static Real Minimum()
      // { return MINDOUBLE; }
      { return 2.2250738585072014e-308; }
//   FREE_CHECK(FloatingPointPrecision)
};

#endif

#endif




#endif
