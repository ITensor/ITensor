//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#ifndef __ITENSOR_QN_H
#define __ITENSOR_QN_H

#include "itensor/global.h"
#include "itensor/arrow.h"
#include "itensor/smallstring.h"
#include "itensor/util/h5/wrap_h5.hpp"

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

//QNum: storage element type for QN
//Represents a number with a Z (integer)
//or Z_M (integer mod M) addition rule
//
// Meaning of mod field:
// mod == 1  => Z addition
// mod >  1  => Z_M addition
// mod <  0  => same as above, fermionic
// mod == 0  => inactive/not used
struct QNum
    {
    using qn_t = int;
    private:
    QNName name_;
    qn_t val_ = 0,
         mod_ = 0;
    public:

    QNum() { }

    explicit 
    QNum(qn_t v) : name_(""), mod_(1) { set(v); }

    QNum(QNName name, qn_t v) : name_(name), mod_(1) { set(v); }

    QNum(qn_t v, qn_t m) : name_(""), mod_(m) { set(v); }

    QNum(QNName name, qn_t v, qn_t m) : name_(name), mod_(m) { set(v); }

    //QNum(std::tuple<std::string,int,int> qv)
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

    QNum
    operator-() const { return QNum(name_,-val_, mod_); }

    explicit operator bool() const { return mod_ != 0; }
    };

//
// QN - quantum number label for IQIndex sectors
//

class QN
    {
    public:
    using qn_t = QNum::qn_t;
    using storage_type = std::array<QNum,QNSize()>;
    using iterator = storage_type::iterator;
    using const_iterator = storage_type::const_iterator;
    private:
    storage_type qvs_{};
    public:

    QN() { }

    explicit
    QN(qn_t q0);

    explicit
    QN(QNum v0,
       QNum v1 = QNum{},
       QNum v2 = QNum{},
       QNum v3 = QNum{}) 
       { 
       addNum(v0);
       if(not v1) return;
       addNum(v1);
       if(not v2) return;
       addNum(v2);
       if(not v3) return;
       addNum(v3);
       }

    explicit operator bool() const { return qvs_.front().mod() != 0; }

    size_t
    size() const { return qvs_.size(); }

    void
    addNum(QNum const& qv);

    QNum 
    num(QNName name) const;

    qn_t //qn_t == int
    val(QNName const& name) const;

    qn_t //qn_t == int
    mod(QNName const& name) const;

    bool
    hasName(QNName const& name) const;

    //1-indexed
    QNum &
    num(size_t n);

    //1-indexed
    QNum const&
    num(size_t n) const;

    QNName
    name(size_t n) const { return num(n).name(); }

    qn_t
    val(size_t n) const { return num(n).val(); }

    qn_t
    mod(size_t n) const { return num(n).mod(); }

    storage_type &
    store() { return qvs_; }

    storage_type const&
    store() const { return qvs_; }

    void
    modAssign(QN const& qo);
    };



/////////////////////////////////
//
// QNum functions
// 
/////////////////////////////////

bool inline
isFermionic(QNum const& qv) { return qv.mod() < 0; }

bool inline
isActive(QNum const& qv) { return qv.mod() != 0; }

QNum&
operator+=(QNum& qva, QNum const& qvb);

QNum&
operator-=(QNum& qva, QNum const& qvb);

QNum&
operator*=(QNum& qva, Arrow dir);

QNum inline
operator+(QNum qva, QNum const& qvb) { return qva += qvb; }

QNum inline
operator-(QNum qva, QNum const& qvb) { return qva -= qvb; }

bool
operator==(QNum const& qva, QNum const& qvb);

bool
operator!=(QNum const& qva, QNum const& qvb);

void
read(std::istream & s, QNum & q);

void
write(std::ostream & s, QNum const& q);


//
// QN functions
// 

bool inline
isFermionic(QN const& q, QNName const& name) 
    { 
    return isFermionic(q.num(name));
    }

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
operator<<(std::ostream & s, QNum const& qv);

inline QNum QN::
num(QNName name) const 
    { 
    for(auto& v : qvs_)
        {
        if(v.name() == name) return v;
        }
    Error("QNum with given name not found");
    return QNum();
    }

inline QNum& QN::
num(size_t n)
    { 
#ifdef DEBUG
    return qvs_.at(n-1); 
#else
    return qvs_[n-1]; 
#endif
    }

inline QNum const& QN::
num(size_t n) const 
    { 
#ifdef DEBUG
    return qvs_.at(n-1); 
#else
    return qvs_[n-1]; 
#endif
    }

//QN inline
//zeroOf(QN q)
//    {
//    for(auto n : range1(QNSize()))
//        {
//        q.num(n).set(0);
//        }
//    return q;
//    }

#ifdef ITENSOR_USE_HDF5
void
h5_write(h5::group parent, std::string const& name, QN const& q);
void
h5_read(h5::group parent, std::string const& name, QN & q);
#endif //ITENSOR_USE_HDF5

} //namespace itensor

#endif
