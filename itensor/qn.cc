#include "itensor/qn.h"
#include "itensor/util/readwrite.h"

namespace itensor {

QN::
QN(qn_t q0)
    {
    qn_[0] = QNVal(q0);
    }

QN::
QN(qn_t q0,
   qn_t q1)
    {
    qn_[0] = QNVal(q0);
    qn_[1] = QNVal(q1);
    }

QN::
QN(qn_t q0,
   qn_t q1,
   qn_t q2)
    {
    qn_[0] = QNVal(q0);
    qn_[1] = QNVal(q1);
    qn_[2] = QNVal(q2);
    }

QN::
QN(qn_t q0,
   qn_t q1,
   qn_t q2,
   qn_t q3)
    {
    qn_[0] = QNVal(q0);
    qn_[1] = QNVal(q1);
    qn_[2] = QNVal(q2);
    qn_[3] = QNVal(q3);
    }

QN::
QN(Args const& args)
    {
    auto hasSz = args.defined("Sz");
    auto start = hasSz ? 1 : 0;

    if(hasSz) qn_[0] = QNVal(args.getInt("Sz"));

    if(args.defined("Nb"))
        {
        qn_[start] = QNVal(args.getInt("Nb"),+1);
        }
    else if(args.defined("Nf"))
        {
        qn_[start] = QNVal(args.getInt("Nf"),-1);
        }
    else if(args.defined("Pf"))
        {
        qn_[start] = QNVal(args.getInt("Pf"),-2);
        }
    }

void QN::
modAssign(QN const& qo)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        qn_[n] = QNVal(qn_[n].val(),qo.qn_[n].mod());
        }
    }

void QNVal::
set(qn_t v)
    { 
    auto m = std::abs(mod_);
    if(m > 1) 
        {
        // Here v is the value we want to compute mod m.
        // If v >= 0 then it's enough to compute v%m.
        // If v < 0 we don't want -(|v|%m) which is what
        // C++'s v%m operation gives. Instead we want
        // to map, say, -1 to m-1 and -m-1 to m-1 etc.
        // Observe that m*|v| > v (since m > 1).
        // So m*|v|+v for negative v will always be > 0
        // and mod'ing this will give the desired behavior.
        val_ = (m*std::abs(v)+v)%m;
        }
    else                   
        {
        val_ = v;
        }
    }


void
operator+=(QNVal& qva, QNVal const& qvb) 
    { 
#ifdef DEBUG
    if(qva.mod() != qvb.mod())
        {
        printfln("qva.mod()=%d, qvb.mod()=%d",qva.mod(),qvb.mod());
        Error("Mismatched mod factors");
        }
#endif
    qva.set(qva.val()+qvb.val());
    }

void
operator-=(QNVal& qva, QNVal const& qvb) 
    { 
    assert(qva.mod() == qvb.mod());
    qva.set(qva.val()-qvb.val());
    }

void
operator*=(QNVal& qva, Arrow dir)
    { 
    qva.set(qva.val() * static_cast<int>(dir));
    }

bool
operator==(QNVal const& qva, QNVal const& qvb)
    {
    assert(qva.mod() == qvb.mod());
    return qva.val() == qvb.val();
    }

bool
operator!=(QNVal const& qva, QNVal const& qvb)
    {
    assert(qva.mod() == qvb.mod());
    return qva.val() != qvb.val();
    }

void
read(std::istream & s, QNVal & q)
    {
    QNVal::qn_t v = 0,
                m = 0;
    itensor::read(s,v);
    itensor::read(s,m);
    q = QNVal(v,m);
    }

void
write(std::ostream & s, QNVal const& q)
    {
    itensor::write(s,q.val());
    itensor::write(s,q.mod());
    }

void
printFull(QN const& q)
    {
    print("QN(");
    for(auto n : range1(QNSize()))
        {
        if(!isActive(q,n)) break;
        if(n > 1) print(",");
        print("{",q(n),",",q.mod(n),"}");
        }
    println(")");
    }

bool
operator==(QN const& qa, QN const& qb)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        if(qa.val0(n).val() != qb.val0(n).val()) return false;
        }
    return true;
    }

bool
operator<(QN const& qa, QN const& qb)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        if(qa.val0(n).val() == qb.val0(n).val()) continue;
        return qa.val0(n).val() < qb.val0(n).val();
        }
    return false;
    }

QN
operator-(QN q)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        q.val0(n) = -q.val0(n);
        }
    return q;
    }

void
operator+=(QN & qa, QN const& qb) 
    { 
    if(!qa) qa.modAssign(qb);
    if(!qb) return;
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        qa.val0(n) += qb.val0(n);
        }
    }

void
operator-=(QN & qa, QN const& qb) 
    { 
    if(!qa) qa.modAssign(qb);
    if(!qb) return;
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        qa.val0(n) -= qb.val0(n);
        }
    }

void
operator*=(QN & qa, Arrow dir)
    { 
    for(size_t n = 0; n < QNSize() && isActive(qa.val0(n)); ++n)
        {
        qa.val0(n) *= dir;
        }
    }

std::ostream& 
operator<<(std::ostream & s, QN const& q)
    {
    if(q.mod(1) == 1 && !isActive(q,2))
        {
        //spin or spinless boson
        s << "QN(" << q(1) << ")";
        }
    else
    if(q.mod(1) == -1 && !isActive(q,2))
        {
        //spinless fermion
        s << "(Nf=" << q(1) << ")";
        }
    else
    if(q.mod(1) == -2 && !isActive(q,2))
        {
        //parity-only spinless fermion
        s << "(Pf=" << q(1) << ")";
        }
    else
    if(q.mod(1) == 1 && q.mod(2) == -1 && !isActive(q,3))
        {
        //electron
        s << "(Sz=" << q(1) << ",Nf=" << q(2) << ")";
        }
    else
    if(q.mod(1) == 1 && q.mod(2) == -2 && !isActive(q,3))
        {
        //"superconducting" electron (parity conservation only)
        s << "(Sz=" << q(1) << ",Pf=" << q(2) << ")";
        }
    else
        {
        //catch-all behavior
        s << "QN(";
        for(auto n : range1(QNSize()))
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
        s << ")";
        }
    return s;
    }

int
paritySign(QN const& q)
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

bool
isFermionic(QN const& q)
    {
    for(size_t n = 0; n < QNSize() && isActive(q.val0(n)); ++n)
        {
        if(isFermionic(q.val0(n))) return true;
        }
    return false;
    }

void
read(std::istream & s, QN & q)
    {
    for(auto& el : q.store())
        itensor::read(s,el);
    }

void
write(std::ostream & s, QN const& q)
    {
    for(auto& el : q.store())
        itensor::write(s,el);
    }
} //namespace itensor
