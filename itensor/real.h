//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_REAL_H
#define __ITENSOR_REAL_H
#include "matrix.h"
#include <limits>
#include "math.h"
#include <cmath>
#include "print.h"
#include "types.h"

#ifndef NAN
#define NAN (std::numeric_limits<Real>::quiet_NaN())
#endif

namespace itensor {

//Real ran1();

static const Real Pi = M_PI;
static const Real Sqrt2 = sqrt(2);
static const Real ISqrt2 = 1.0/sqrt(2);
//static const Real Sqrt3 = sqrt(3);
//static const Real ISqrt3 = 1.0/sqrt(3);

template <typename T>
T sqr(T x) { return x*x; }

static const Real ApproxReal_Accuracy = 1E-12;

struct ApproxReal
    {
    Real r;

    //Default constructed to NAN 
    //to signal initialization errors
    ApproxReal() : r(NAN) {}

    ApproxReal(Real _r) : r(_r) {}

    bool 
    operator==(const ApproxReal& other) const
        { return fabs(r-other.r) <= ApproxReal_Accuracy; }

    bool 
    operator!=(const ApproxReal& other) const
        { return fabs(r-other.r) > ApproxReal_Accuracy; }

    bool
    operator<(const ApproxReal& other) const
        { return other.r-r > ApproxReal_Accuracy; }

    ApproxReal& 
    operator+=(const ApproxReal &A)
        { r += A.r; return *this; }

    ApproxReal& 
    operator+=(Real a)
        { r += a; return *this; }

    };

static const Real maxlogdouble = log(std::numeric_limits<double>::max());

static const Real LogNumber_Accuracy = 1E-12;

class TooBigForReal : public ITError
    {
public:
    typedef ITError
    Parent;

    TooBigForReal(const std::string& message) 
        : Parent(message)
        { }
    };

class TooSmallForReal : public ITError
    {
public:
    typedef ITError
    Parent;

    TooSmallForReal(const std::string& message) 
        : Parent(message)
        { }
    };

//
// LogNumber
//

//Stores a real number r as lognum_ = log(|r|) and sign_ = sgn(r)
class LogNumber
    {
    public:

    Real 
    logNum() const { return lognum_; }

    int 
    sign() const { return sign_; }

    bool
    isZero() const { return sign_ == 0; }

    bool 
    isRealZero() const
        { return (sign_ == 0 || lognum_ < -maxlogdouble); }

    bool 
    isFiniteReal() const 
        { return (lognum_ < maxlogdouble && lognum_ > -maxlogdouble); }

    bool 
    isTooBigForReal() const 
        { return (lognum_ > maxlogdouble); }

    bool 
    isTooSmallForReal() const 
        { return (lognum_ < -maxlogdouble); }

    bool friend inline
    isnan(const LogNumber& L) { return std::isnan(L.lognum_); }

    //Default constructed to NAN 
    //to catch initialization errors
    LogNumber() 
        : 
        lognum_(NAN), 
        sign_(1) 
        { }

    explicit
    LogNumber(Real r)
        {
        if(r == 0)
            { 
            sign_ = 0;  
            lognum_ = 0; 
            }
        else 
        if(r < 0)
            { 
            sign_ = -1; 
            lognum_ = log(-r); 
            }
        else
            { 
            sign_ = 1;  
            lognum_ = log(r); 
            }
        }

    LogNumber(Real lognum, int sign) 
        : 
        lognum_(lognum), 
        sign_(sign) 
        { 
#ifdef DEBUG
        if(!(sign == -1 || sign == 0 || sign == 1))
            Error("sign should be -1, 0, or 1");
#endif
        } 

    void 
    read(std::istream& s) 
        { 
        s.read((char*) &lognum_,sizeof(lognum_)); 
        s.read((char*) &sign_,sizeof(sign_)); 
        }
    void 
    write(std::ostream& s) const 
        { 
        s.write((char*) &lognum_,sizeof(lognum_)); 
        s.write((char*) &sign_,sizeof(sign_)); 
        }

    //operator Real() const	too easy to misuse accidentally
    Real 
    real() const
        {
        if(sign_ == 0) return 0;
#ifdef DEBUG
        if(lognum_ > maxlogdouble)
            { 
            println("lognum_ = ",lognum_);
            throw TooBigForReal("LogNumber too big to convert to Real");
            }
        if(lognum_ < -maxlogdouble)
            { 
            println("lognum_ = ",lognum_);
            throw TooSmallForReal("LogNumber too small to convert to Real");
            }
#endif
        return sign_ * exp(lognum_);
        }

    Real
    real0() const
        {
        if(isRealZero()) return 0.0;
        return real();
        }

    LogNumber&
    operator+=(const LogNumber& other)
        {
        try {
            *this = 
                LogNumber(this->real()+other.real());
            }
        catch(const ITError& e)
            {
            Error("Could not convert to real in LogNumber::operator+=");
            }
        return *this;
        }


    bool 
    operator==(const LogNumber& other) const
        { return (sign_ == other.sign_) && (lognum_ == other.lognum_); }

    bool 
    operator!=(const LogNumber& other) const
        { return (sign_ != other.sign_) || (lognum_ != other.lognum_); }

    bool 
    approxEquals(const LogNumber& other) const
        { return (sign_ == other.sign_) && (fabs(lognum_-other.lognum_) < LogNumber_Accuracy); }

    void
    negate() { sign_ *= -1; }

    LogNumber& 
    operator*=(const LogNumber& other)
        {
        sign_ *= other.sign_;
        lognum_ += other.lognum_;
        return *this;
        }

    LogNumber& 
    operator*=(Real other)
        {
        if(other == -1)	// special case for efficiency
            {
            sign_ *= -1;
            return *this;
            }
        return *this *= LogNumber(other);
        }

    LogNumber& 
    operator/=(const LogNumber& other)
        {
#ifdef DEBUG
        if(other.sign_ == 0) Error("divide by zero in LogNumber");
#endif
        sign_ *= other.sign_;
        lognum_ -= other.lognum_;
        return *this;
        }

    LogNumber& 
    operator/=(Real other) { return *this /= LogNumber(other); }

    LogNumber 
    operator-() const 
        { LogNumber res(*this); res.sign_ *= -1; return res; }

    LogNumber 
    operator/(const LogNumber& other) const
        { LogNumber res(*this); res /= other; return res; }
    LogNumber 
    operator*(const LogNumber& other) const
        { LogNumber res(*this); res *= other; return res; }

    bool 
    operator<(const LogNumber& other) const
        {
        if(sign_ != other.sign_)
            { return sign_ < other.sign_; }
        else if(sign_ > 0)  //this and other are positive
            { return lognum_ < other.lognum_; }
        else if(sign_ == 0) //this and other are zero
            { return false; }
        //this and other are negative
        return lognum_ > other.lognum_;
        }

    bool 
    operator<=(const LogNumber& other) const
        {
        return (operator<(other) || operator==(other));
        }

    bool 
    operator>=(const LogNumber& other) const
        {
        return !(operator<(other));
        }

    bool 
    operator>(const LogNumber& other) const
        {
        return (!operator<(other)) && (!operator==(other));
        }

    void
    swap(LogNumber& other)
        {
        Real sl = lognum_;
        lognum_ = other.lognum_;
        other.lognum_ = sl;

        int si = sign_;
        sign_ = other.sign_;
        other.sign_ = si;
        }

    bool 
    magnitudeLessThan(const LogNumber& other) const
        {
        if(sign_ == 0) return other.sign_ != 0;
        if(other.sign_ == 0) return false;
        return lognum_ < other.lognum_;
        }

    const LogNumber&
    pow(Real p)
        {
        lognum_ *= p;
        return *this;
        }


    private:

    /////////////
    //
    // Data Members
    //

    //Log of the magnitude of 
    //the number represented.
    Real lognum_;

    int sign_;

    //
    /////////////

    };

LogNumber inline
sqrt(LogNumber L)
    {
    if(L.sign() < 0) 
        Error("Negative LogNumber in sqrt");
    return L.pow(0.5);
    }

inline 
std::ostream& 
operator<<(std::ostream& s, const LogNumber& N)
    {
    s << "LogNumber(" << N.logNum() << ",";
    if(N.sign() == 0) s << "0)";
    else           s << (N.sign() > 0 ? "+)" : "-)");
    return s;
    }

}; //namespace itensor

#endif
