#ifndef __ITENSOR_REAL_H
#define __ITENSOR_REAL_H
#include "matrix.h"
#include <limits>
#include "types.h"

static const Real Pi = M_PI;
static const Real Sqrt2 = sqrt(2);
static const Real ISqrt2 = 1.0/sqrt(2);
static const Real Sqrt3 = sqrt(3);
static const Real ISqrt3 = 1.0/sqrt(3);

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

static Real maxlogdouble = log(std::numeric_limits<double>::max());

static const Real LogNumber_Accuracy = 1E-12;

//Stores a real number r as lognum_ = log(|r|) and sign_ = sgn(r)
class LogNumber
{
    Real lognum_;
    int sign_;
public:

    Real logNum() const { return lognum_; }
    int sign() const { return sign_; }
    inline bool isRealZero() const
        { return (sign_ == 0 || lognum_ < -maxlogdouble); }
    inline bool isOne() const { return (lognum_ == 0 && sign_ == 1); }
    inline bool isFinite() const 
        { return (lognum_ < maxlogdouble && lognum_ > -maxlogdouble); }
    inline bool isNan() const 
        { return lognum_ != lognum_; }

    //Default is Real(LogNum()) == 1
    LogNumber() : lognum_(0), sign_(1) { }

    LogNumber(Real r)
	{
        if(r == 0)
            { sign_ = 0;  lognum_ = 0; }
        else if(r > 0)
            { sign_ = 1;  lognum_ = log(r); }
        else
            { sign_ = -1; lognum_ = log(-r); }
	}

    LogNumber(Real lognum, int sign) : lognum_(lognum), sign_(sign) { } 

    inline void read(std::istream& s) 
        { s.read((char*)this,sizeof(this)); }
    inline void write(std::ostream& s) const 
        { s.write((char*)this,sizeof(this)); }

    //operator Real() const	too easy to misuse accidentally
    Real real() const
	{
        if(sign_ == 0) return 0;
#ifndef DNDEBUG
        if(lognum_ > maxlogdouble)
        { 
            Print(lognum_);
            Error("LogNumber too big to convert to Real"); 
        }
        if(lognum_ < -maxlogdouble)
        { 
            Print(lognum_);
            Error("LogNumber too small to convert to Real"); 
        }
#endif
        return sign_ * exp(lognum_);
	}

    inline bool operator==(const LogNumber& other) const
    { return (sign_ == other.sign_) && (lognum_ == other.lognum_); }

    inline bool approxEquals(const LogNumber& other) const
    { return (sign_ == other.sign_) && (fabs(lognum_-other.lognum_) < LogNumber_Accuracy); }

    LogNumber& operator*=(const LogNumber& other)
	{
        sign_ *= other.sign_;
        lognum_ += other.lognum_;
        return *this;
	}

    LogNumber& operator*=(Real other)
	{
        if(other == -1)	// special case for efficiency
        {
            sign_ *= -1;
            return *this;
        }
        return *this *= LogNumber(other);
	}

    LogNumber& operator/=(const LogNumber& other)
	{
        DO_IF_DEBUG(if(other.sign_ == 0) Error("divide by zero in LogNumber");)
        sign_ *= other.sign_;
        lognum_ -= other.lognum_;
        assert(lognum_ < maxlogdouble);
        assert(lognum_ > -maxlogdouble);
        return *this;
	}

    LogNumber& operator/=(Real other) { return *this /= LogNumber(other); }

    LogNumber operator/(const LogNumber& other)
	{ LogNumber res(*this); res /= other; return res; }
    LogNumber operator*(const LogNumber& other)
	{ LogNumber res(*this); res *= other; return res; }

    bool operator<(const LogNumber& other) const
	{
        if(sign_ != other.sign_)
            { return sign_ < other.sign_; }
        else if(sign_ == 0)
            { return false; }
        else if(sign_ > 0)
            { return lognum_ < other.lognum_; }
        return lognum_ > other.lognum_;
	}

    bool magnitudeLessThan(const LogNumber& other) const
	{
        if(sign_ == 0) return other.sign_ != 0;
        if(other.sign_ == 0) return false;
        return lognum_ < other.lognum_;
	}

    friend inline ostream& operator<<(ostream& s, const LogNumber& N)
    {
        s << "LogNumber(" << N.logNum() << ",";
        if(N.sign() == 0) s << "0)";
        else           s << (N.sign() > 0 ? "+)" : "-)");
        return s;
    }

};

#endif
