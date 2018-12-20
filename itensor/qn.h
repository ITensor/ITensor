//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QN_H
#define __ITENSOR_QN_H

#include "itensor/global.h"
#include "itensor/arrow.h"

namespace itensor {

//
// QN
//
// Quantum number label for IQIndexes.
//

class QN;

using QNName = SmallString;

size_t inline constexpr
QNSize() { return 4ul; }

//QNVal: storage element type for QN
//Represents a number with a Z (integer)
//or Z_M (integer mod M) addition rule
//
// Meaning of mod field:
// mod == 1  => Z addition
// mod >  1  => Z_M addition
// mod <  0  => same as above, fermionic
// mod == 0  => inactive/not used
struct QNVal
    {
    using qn_t = int;
    private:
    QNName name_;
    qn_t val_ = 0,
         mod_ = 0;
    public:

    QNVal() { }

    explicit QNVal(qn_t v) : name_(""), val_(v), mod_(1) { }

    QNVal(QNName name, qn_t v) : name_(name), val_(v), mod_(1) { }

    QNVal(QNName name, qn_t v, qn_t m) : name_(name), mod_(m) { set(v); }

    QNVal(std::initializer_list<qn_t> qv)
        {
        if(qv.size() != 3) Error("initializer_list arg to QNVal must have three elements");
        name_ = *(qv.begin());
        mod_ = *(qv.begin()+2);
        set(*(qv.begin()+1));
        }

    QNName
    name() const { return name_; }

    qn_t
    mod() const { return mod_; }

    qn_t
    val() const { return val_; }

    void
    set(qn_t v);

    QNVal
    operator-() const { return QNVal(name_,-val_, mod_); }
    };

//
// QN - quantum number label for IQIndex sectors
//

class QN
    {
    public:
    using qn_t = QNVal::qn_t;
    using storage_type = std::array<QNVal,QNSize()>;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;
    private:
    storage_type qn_{};
    public:

    QN() { }

    //// Takes named Args:
    //// QN({"Sz=",-1,"Nf=",2})
    //explicit
    //QN(Args const& args);

    explicit
    QN(qn_t q0);

    QN(qn_t q0,
       qn_t q1);

    QN(qn_t q0,
       qn_t q1,
       qn_t q2);

    QN(qn_t q0,
       qn_t q1,
       qn_t q2,
       qn_t q3);

    template <typename T, typename... Rest>
    QN(const char* name1,
       T const& t1, 
       Rest const&... rest)
      : QN(Args(name1,t1,rest...))
        { }

    // Takes QNVal arguments,
    // specifying both a qn value in each
    // sector and a mod factor
    // Can call as as QN({0,2},{1,2})
    explicit
    QN(QNVal v0,
       QNVal v1 = QNVal{},
       QNVal v2 = QNVal{},
       QNVal v3 = QNVal{}) 
     : qn_{{v0,v1,v2,v3}}  //TODO: make sure this is sorted
       { }

    explicit
    QN(std::initializer_list<qn_t> qv)
        {
        qn_[0] = qv;
        }

    explicit operator bool() const { return qn_.front().mod() != 0; }

    size_t
    size() const { return qn_.size(); }

//    qn_t
//    operator[](size_t n) const
//        { 
//#ifdef DEBUG
//        return qn_.at(n).val(); 
//#else
//        return qn_[n].val(); 
//#endif
//        }
//
//    //1-indexed
//    qn_t
//    operator()(size_t n) const { return operator[](n-1); }
//
//    //1-indexed
//    qn_t
//    mod(size_t n) const { return qn_.at(n-1).mod(); }

    //0-indexed
    QNVal &
    val0(size_t n)
        { 
#ifdef DEBUG
        return qn_.at(n); 
#else
        return qn_[n]; 
#endif
        }

    //0-indexed
    QNVal const&
    val0(size_t n) const 
        { 
#ifdef DEBUG
        return qn_.at(n); 
#else
        return qn_[n]; 
#endif
        }

    void
    modAssign(QN const& qo);

    storage_type &
    store() { return qn_; }

    storage_type const&
    store() const { return qn_; }
    };


//
// QNVal functions
// 

bool inline
isFermionic(QNVal const& qv) { return qv.mod() < 0; }

bool inline
isActive(QNVal const& qv) { return qv.mod() != 0; }

void
operator+=(QNVal& qva, QNVal const& qvb);

void
operator-=(QNVal& qva, QNVal const& qvb);

void
operator*=(QNVal& qva, Arrow dir);

bool
operator==(QNVal const& qva, QNVal const& qvb);

bool
operator!=(QNVal const& qva, QNVal const& qvb);

void
read(std::istream & s, QNVal & q);

void
write(std::ostream & s, QNVal const& q);


//
// QN functions
// 

bool
operator==(QN const& qa, QN const& qb);

bool inline
operator!=(QN const& qa, QN const& qb) { return !operator==(qa,qb); }

bool
operator<(QN const& qa, QN const& qb);

QN
operator-(QN q);

void
operator+=(QN & qa, QN const& qb);

void
operator-=(QN & qa, QN const& qb);

void
operator*=(QN & qa, Arrow dir);

QN inline
operator+(QN qa, QN const& qb) { qa += qb; return qa; }

QN inline
operator-(QN qa, QN const& qb) { qa -= qb; return qa; }

QN inline
operator*(QN q, Arrow dir) { q *= dir; return q; }

QN inline
operator*(Arrow dir, QN q) { q *= dir; return q; }

std::ostream& 
operator<<(std::ostream & s, QN const& q);

void
read(std::istream & s, QN & q);

void
write(std::ostream & s, QN const& q);

void
printFull(QN const& q);

} //namespace itensor

#endif
