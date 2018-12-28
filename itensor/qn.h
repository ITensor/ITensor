//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QN_H
#define __ITENSOR_QN_H

#include "itensor/global.h"
#include "itensor/arrow.h"
#include "itensor/smallstring.h"

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

    explicit 
    QNVal(qn_t v) : name_(""), mod_(1) { set(v); }

    QNVal(QNName name, qn_t v) : name_(name), mod_(1) { set(v); }

    QNVal(qn_t v, qn_t m) : name_(""), mod_(m) { set(v); }

    QNVal(QNName name, qn_t v, qn_t m) : name_(name), mod_(m) { set(v); }

    //QNVal(std::tuple<std::string,int,int> qv)
    //    {
    //    name_ = QNName(std::get<0>(qv));
    //    mod_ = std::get<2>(qv);
    //    auto value = std::get<1>(qv);
    //    set(value);
    //    }

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

    explicit operator bool() const { return mod_ != 0; }
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
    storage_type qvs_{};
    public:

    QN() { }

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

    explicit
    QN(QNVal v0,
       QNVal v1 = QNVal{},
       QNVal v2 = QNVal{},
       QNVal v3 = QNVal{}) 
       { 
       addVal(v0);
       if(not v1) return;
       addVal(v1);
       if(not v2) return;
       addVal(v2);
       if(not v3) return;
       addVal(v3);
       }

    explicit operator bool() const { return qvs_.front().mod() != 0; }

    size_t
    size() const { return qvs_.size(); }

    void
    addVal(QNVal const& qv);

    qn_t //qn_t == int
    getVal(QNName const& name) const;

    storage_type &
    store() { return qvs_; }

    storage_type const&
    store() const { return qvs_; }

    //0-indexed
    QNVal &
    val0(size_t n);

    //0-indexed
    QNVal const&
    val0(size_t n) const;

    void
    modAssign(QN const& qo);
    };

/////////////////////////////////
//
// QNVal functions
// 
/////////////////////////////////

bool inline
isFermionic(QNVal const& qv) { return qv.mod() < 0; }

bool inline
isActive(QNVal const& qv) { return qv.mod() != 0; }

QNVal&
operator+=(QNVal& qva, QNVal const& qvb);

QNVal&
operator-=(QNVal& qva, QNVal const& qvb);

QNVal&
operator*=(QNVal& qva, Arrow dir);

QNVal inline
operator+(QNVal qva, QNVal const& qvb) { return qva += qvb; }

QNVal inline
operator-(QNVal qva, QNVal const& qvb) { return qva -= qvb; }

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
operator==(QN qa, QN const& qb);

bool inline
operator!=(QN const& qa, QN const& qb) { return !operator==(qa,qb); }

bool
operator<(QN const& qa, QN const& qb);

QN
operator-(QN q);

QN&
operator+=(QN & qa, QN const& qb);

QN&
operator-=(QN & qa, QN const& qb);

QN&
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

std::ostream& 
operator<<(std::ostream & s, QNVal const& qv);

inline QNVal& QN::
val0(size_t n)
    { 
#ifdef DEBUG
    return qvs_.at(n); 
#else
    return qvs_[n]; 
#endif
    }

inline QNVal const& QN::
val0(size_t n) const 
    { 
#ifdef DEBUG
    return qvs_.at(n); 
#else
    return qvs_[n]; 
#endif
    }

} //namespace itensor

#endif
