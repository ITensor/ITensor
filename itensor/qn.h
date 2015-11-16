//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_QN_H
#define __ITENSOR_QN_H

namespace itensor {

//
// QN
//
// Quantum number label for IQIndexes.
//
//

int
qnrem(int v, int m);

struct QNVal
    {
    using qn_t = int;
    private:
    qn_t val_ = 0,
         mod_ = 0;
    public:

    QNVal() { }

    QNVal(qn_t v) : val_(v), mod_(1) { }

    QNVal(qn_t v, qn_t m) : mod_(m) { set(v); }

    qn_t
    mod() const { return mod_; }

    qn_t
    val() const { return val_; }

    void
    set(qn_t v)
        { 
        if(std::abs(mod_) > 1) 
            {
            val_ = std::abs(v%mod_);
            }
        else                   
            {
            val_ = v;
            }
        }

    QNVal&
    operator-() { val_ = -val_; return *this; }
    };

size_t inline constexpr
QNSize() { return 4ul; }

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

    template<typename... VArgs>
    explicit
    QN(QNVal v0,
        VArgs&&... vals)
      : qn_{{v0,QNVal(vals)...}}
        { 
        static_assert(1+sizeof...(VArgs) <= QNSize(),"Too many arguments to QN constructor");
        }

    template<typename... Qs>
    explicit
    QN(qn_t q0,
        Qs&&... qs)
      : qn_{{QNVal(q0),QNVal(qs)...}}
        { 
        static_assert(1+sizeof...(Qs) <= QNSize(),"Too many arguments to QN constructor");
        }

    explicit operator bool() const { return qn_.front().mod() != 0; }

//    qn_t &
//    operator[](size_t n) 
//        { 
//#ifdef DEBUG
//        return qn_.at(n).val(); 
//#else
//        return qn_[n].val(); 
//#endif
//        }

    qn_t
    operator[](size_t n) const
        { 
#ifdef DEBUG
        return qn_.at(n).val(); 
#else
        return qn_[n].val(); 
#endif
        }

//    //1-indexed
//    qn_t & 
//    operator()(size_t n) { return operator[](n-1); }

    //1-indexed
    qn_t
    operator()(size_t n) const { return operator[](n-1); }

    //1-indexed
    qn_t
    mod(size_t n) const { return qn_.at(n-1).mod(); }

    size_t
    size() const { return qn_.size(); }

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
    modAssign(QN const& qo)
        {
        for(size_t n = 0; n < QNSize(); ++n)
            {
            qn_[n] = QNVal(qn_[n].val(),qo.qn_[n].mod());
            }
        }


    storage_type &
    store() { return qn_; }

    storage_type const&
    store() const { return qn_; }
    };

bool inline
isFermionic(QNVal const& qv) { return qv.mod() < 0; }

bool inline
isActive(QNVal const& qv) { return qv.mod() != 0; }

void inline
operator+=(QNVal& qva, QNVal const& qvb) 
    { 
    assert(qva.mod() == qvb.mod());
    qva.set(qva.val()+qvb.val());
    }

void inline
operator-=(QNVal& qva, QNVal const& qvb) 
    { 
    assert(qva.mod() == qvb.mod());
    qva.set(qva.val()-qvb.val());
    }

void inline
operator*=(QNVal& qva, Arrow dir)
    { 
    qva.set(qva.val() * static_cast<int>(dir));
    }

bool inline
operator==(QNVal const& qva, QNVal const& qvb)
    {
    assert(qva.mod() == qvb.mod());
    return qva.val() == qvb.val();
    }

bool inline
operator!=(QNVal const& qva, QNVal const& qvb)
    {
    assert(qva.mod() == qvb.mod());
    return qva.val() != qvb.val();
    }

bool
isActive(QN const& q, size_t n) { return isActive(q.val0(n-1)); }

bool inline
isFermionic(QN const& q, size_t n) { return q.mod(n) < 0; }

bool inline
operator==(QN const& qa, QN const& qb)
    {
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        if(qa.val0(n).val() != qb.val0(n).val()) return false;
        }
    return true;
    }

bool inline
operator!=(QN const& qa, QN const& qb) { return !operator==(qa,qb); }

void inline
operator+=(QN & qa, QN const& qb) 
    { 
    if(!qa) qa.modAssign(qb); 
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        qa.val0(n) += qb.val0(n);
        }
    }

void inline
operator-=(QN & qa, QN const& qb) 
    { 
    if(!qa) qa.modAssign(qb);
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        qa.val0(n) -= qb.val0(n);
        }
    }

void inline
operator*=(QN & qa, Arrow dir)
    { 
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        qa.val0(n) *= dir;
        }
    }

QN inline
operator+(QN qa, QN const& qb) { qa += qb; return qa; }
QN inline
operator-(QN qa, QN const& qb) { qa -= qb; return qa; }
QN inline
operator*(QN q, Arrow dir) { q *= dir; return q; }
QN inline
operator*(Arrow dir, QN q) { q *= dir; return q; }

inline std::ostream& 
operator<<(std::ostream & s, QN const& q)
    {
    s << "QN(";
    if(q.mod(1) == 1 && !isActive(q,2))
        {
        //spin or spinless boson
        s << q(1);
        }
    else
    if(q.mod(1) == -1 && !isActive(q,2))
        {
        //spinless fermion
        s << "Nf=" << q(1);
        }
    else
    if(q.mod(1) == 1 && q.mod(2) == -1 && !isActive(q,3))
        {
        //electron
        s << "Sz=" << q(1) << ",Nf=" << q(2);
        }
    else
    if(q.mod(1) == 1 && q.mod(2) == -2 && !isActive(q,3))
        {
        //"superconducting" electron (parity conservation only)
        s << "Sz=" << q(1) << ",Pf=" << q(2);
        }
    else
        {
        //catch-all behavior
        for(auto n : count1(QNSize()))
            {
            if(!isActive(q,n)) break;
            if(n > 1) s << ",";
            if(q.mod(n) != 1)
                {
                s << "{" << q(n) << "," << q.mod(n) << "}";
                }
            else
                {
                s << q(n);
                }
            }
        }
    return s << ")";
    }

template<typename T> 
int 
sgn(T val) 
    {
    return (T(0) < val) - (val < T(0));
    }

//Sz in units of spin 1/2
QN inline
spin(int Sz) { return QN(Sz); }

QN inline
boson(int Nb) { return QN(Nb); }

//Sz in units of spin 1/2
QN inline
spinboson(int Sz, int Nb) { return QN(Sz,Nb); }

QN inline
fermion(int Nf) { return QN(QNVal(Nf,-1)); }

//Sz in units of spin 1/2
QN inline
electron(int Sz, int Nf) { return QN(QNVal(Sz),QNVal(Nf,-1)); }

//QN conserving electron spin and parity, not total charge
//Sz in units of spin 1/2
QN inline
elparity(int Sz, int Pf) { return QN(QNVal(Sz),QNVal(Pf,-2)); }

int inline
parity(QN const& q)
    {
    int p = 1;
    for(size_t n = 0; n < QNSize() && isActive(q.val0(n)); ++n)
        {
        if(isFermionic(q.val0(n)) && std::abs(q[n])%2==1)
            {
            p *= -1;
            }
        }
    return p;
    }

bool inline
isFermionic(QN const q) { return parity(q) == -1; }

void inline
read(std::istream & s, QNVal & q)
    {
    QNVal::qn_t v = 0,
                m = 0;
    itensor::read(s,v);
    itensor::read(s,m);
    q = QNVal(v,m);
    }

void inline
write(std::ostream & s, QNVal const& q)
    {
    itensor::write(s,q.val());
    itensor::write(s,q.mod());
    }

void inline
read(std::istream & s, QN & q)
    {
    for(auto& el : q.store())
        itensor::read(s,el);
    }

void inline
write(std::ostream & s, QN const& q)
    {
    for(auto& el : q.store())
        itensor::write(s,el);
    }

} //namespace itensor

#undef DEF_NMAX

#endif
