//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_REAL_H
#define __ITENSOR_REAL_H
#include <cmath>
#include "math.h"
#include "itensor/types.h"
#include "itensor/util/print.h"
#include "itensor/util/error.h"

namespace itensor {

const Real Pi = 3.14159265358979323846;
const Real Sqrt2 = std::sqrt(2.);
const Real ISqrt2 = 1./std::sqrt(2.);

template <typename T>
T 
sqr(T x) { return x*x; }

const Real maxlogdouble = log(std::numeric_limits<double>::max());

const Real LogNum_Accuracy = 1E-12;

class TooBigForReal : public ITError
    {
    public:
    using Parent = ITError;

    TooBigForReal(const std::string& message) 
        : Parent(message)
        { }
    };

class TooSmallForReal : public ITError
    {
    public:
    using Parent = ITError;

    TooSmallForReal(const std::string& message) 
        : Parent(message)
        { }
    };

bool inline
equal(Real x, Real y, Real eps = 1E-12)
    {
    auto ax = std::fabs(x);
    auto ay = std::fabs(y);
    auto scale = (ax < ay ? ay : ax);
    return std::fabs(x-y) <= scale*eps;
    }

//
// LogNum
//

//Stores a real number r as lognum_ = log(|r|) and sign_ = sgn(r)
class LogNum
    {
    //Log of the magnitude of 
    //the number represented.
    Real lognum_;
    int sign_;
    public:

    //Default constructed to NAN 
    //to catch initialization errors
    LogNum() 
        : 
        lognum_(NAN), 
        sign_(1) 
        { }

    explicit
    LogNum(Real r)
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

    LogNum(Real lognum, int sign) 
        : 
        lognum_(lognum), 
        sign_(sign) 
        { 
#ifdef DEBUG
        if(!(sign == -1 || sign == 0 || sign == 1))
            Error("sign should be -1, 0, or 1");
#endif
        } 

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
    isnan(const LogNum& L) { return std::isnan(L.lognum_); }


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
            throw TooBigForReal("LogNum too big to convert to Real");
            }
        if(lognum_ < -maxlogdouble)
            { 
            println("lognum_ = ",lognum_);
            throw TooSmallForReal("LogNum too small to convert to Real");
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

    LogNum&
    operator+=(const LogNum& other)
        {
        try {
            *this = 
                LogNum(this->real()+other.real());
            }
        catch(const ITError& e)
            {
            Error("Could not convert to real in LogNum::operator+=");
            }
        return *this;
        }


    bool 
    operator==(const LogNum& other) const
        { return (sign_ == other.sign_) && (lognum_ == other.lognum_); }

    bool 
    operator!=(const LogNum& other) const
        { return (sign_ != other.sign_) || (lognum_ != other.lognum_); }

    bool 
    approxEquals(const LogNum& other) const
        { return (sign_ == other.sign_) && (std::fabs(lognum_-other.lognum_) < LogNum_Accuracy); }

    void
    negate() { sign_ *= -1; }

    LogNum& 
    operator*=(const LogNum& other)
        {
        sign_ *= other.sign_;
        lognum_ += other.lognum_;
        return *this;
        }

    LogNum& 
    operator*=(Real other)
        {
        if(other == -1)	// special case for efficiency
            {
            sign_ *= -1;
            return *this;
            }
        return *this *= LogNum(other);
        }

    LogNum& 
    operator/=(const LogNum& other)
        {
#ifdef DEBUG
        if(other.sign_ == 0) Error("divide by zero in LogNum");
#endif
        sign_ *= other.sign_;
        lognum_ -= other.lognum_;
        return *this;
        }

    LogNum& 
    operator/=(Real other) { return *this /= LogNum(other); }

    LogNum 
    operator-() const 
        { LogNum res(*this); res.sign_ *= -1; return res; }

    LogNum 
    operator/(const LogNum& other) const
        { LogNum res(*this); res /= other; return res; }
    LogNum 
    operator*(const LogNum& other) const
        { LogNum res(*this); res *= other; return res; }

    bool 
    operator<(const LogNum& other) const
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
    operator<=(const LogNum& other) const
        {
        return (operator<(other) || operator==(other));
        }

    bool 
    operator>=(const LogNum& other) const
        {
        return !(operator<(other));
        }

    bool 
    operator>(const LogNum& other) const
        {
        return (!operator<(other)) && (!operator==(other));
        }

    void
    swap(LogNum& other)
        {
        Real sl = lognum_;
        lognum_ = other.lognum_;
        other.lognum_ = sl;

        int si = sign_;
        sign_ = other.sign_;
        other.sign_ = si;
        }

    bool 
    magnitudeLessThan(const LogNum& other) const
        {
        if(sign_ == 0) return other.sign_ != 0;
        if(other.sign_ == 0) return false;
        return lognum_ < other.lognum_;
        }

    const LogNum&
    pow(Real p)
        {
        lognum_ *= p;
        return *this;
        }
    };

//For backwards compatibility:
using LogNumber = LogNum;

LogNum inline
operator*(LogNum L, Real r) { L *= r; return L; }
LogNum inline
operator*(Real r, LogNum L) { L *= r; return L; }

LogNum inline
sqrt(LogNum L)
    {
    if(L.sign() < 0) 
        Error("Negative LogNum in sqrt");
    return L.pow(0.5);
    }

inline std::ostream& 
operator<<(std::ostream& s, LogNum const& N)
    {
    s << "LogNum(" << N.logNum() << ",";
    if(N.sign() == 0) s << "0)";
    else              s << (N.sign() > 0 ? "+)" : "-)");
    return s;
    }

} //namespace itensor

#endif
