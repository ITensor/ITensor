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
    add(QNVal(q0))
    add(QNVal(q1))
    }

QN::
QN(qn_t q0,
   qn_t q1,
   qn_t q2)
    {
    add(QNVal(q0))
    add(QNVal(q1))
    add(QNVal(q2))
    }

QN::
QN(qn_t q0,
   qn_t q1,
   qn_t q2,
   qn_t q3)
    {
    add(QNVal(q0))
    add(QNVal(q1))
    add(QNVal(q2))
    add(QNVal(q3))
    }

//QN::
//QN(Args const& args)
//    {
//    }

void QN::
modAssign(QN const& qo)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        qn_[n] = QNVal(qn_[n].name(),qn_[n].val(),qo.qn_[n].mod());
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
checkCompatible(QNVal const& qva, QNVal const& qvb) 
    {
#ifdef DEBUG
    if(qva.name() != qvb.name())
        {
        printfln("qva.name()=%s, qvb.name()=%s",qva.name(),qvb.name());
        Error("Mismatched QNVal names");
        }
    if(qva.mod() != qvb.mod())
        {
        printfln("qva.mod()=%d, qvb.mod()=%d",qva.mod(),qvb.mod());
        Error("Mismatched mod factors");
        }
#endif
    }


void
operator+=(QNVal& qva, QNVal const& qvb) 
    { 
    checkCompatible(qva,qvb);
    qva.set(qva.val()+qvb.val());
    }

void
operator-=(QNVal& qva, QNVal const& qvb) 
    { 
    checkCompatible(qva,qvb);
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
    checkCompatible(qva,qvb);
    return qva.val() == qvb.val();
    }

bool
operator!=(QNVal const& qva, QNVal const& qvb)
    {
    checkCompatible(qva,qvb);
    return qva.val() != qvb.val();
    }

void
read(std::istream & s, QNVal & q)
    {
    QNName n;
    QNVal::qn_t v = 0,
                m = 0;
    itensor::read(s,n);
    itensor::read(s,v);
    itensor::read(s,m);
    q = QNVal(n,v,m);
    }

void
write(std::ostream & s, QNVal const& q)
    {
    itensor::write(s,q.name());
    itensor::write(s,q.val());
    itensor::write(s,q.mod());
    }

void
printFull(QN const& q)
    {
    print("QN(");
    for(auto& v : q.store())
        {
        if(!isActive(v)) break;
        if(n > 0) print(",");
        print("{\"",v.name(),"\",",v.val(),",",v.mod(),"}");
        }
    println(")");
    }

bool
operator==(QN const& qa, QN const& qb)
    {
    for(auto n : range(QNSize()))
        {
        if(qa.val0(n).val() != qb.val0(n).val()) return false;
        }
    return true;
    }

bool
operator<(QN const& qa, QN const& qb)
    {
    //TODO:
    for(auto n : range(QNSize()))
        {
        if(qa.val0(n).val() == qb.val0(n).val()) continue;
        return qa.val0(n).val() < qb.val0(n).val();
        }
    return false;
    }

QN
operator-(QN q)
    {
    for(auto& v : q.store())
        {
        v = -v;
        }
    return q;
    }

void
operator+=(QN & qa, QN const& qb) 
    { 
    //TODO:
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
    //TODO:
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
    for(auto& v : qa.store())
        {
        v *= dir;
        }
    }

std::ostream& 
operator<<(std::ostream & s, QN const& q)
    {
    s << "QN(";
    for(auto& v : q.store())
        {
        if(!isActive(v)) break;
        if(v.name() != "")
            {
            s << "{\"" << v.name() << "\"," << v.val();
            }
        else if(v.mod() != 1)
            {
            s << "{" << v.val();
            }
        else
            {
            s << v.val();
            }
        if(v.mod() != 1)
            {
            s << "," << v.mod() << "}";
            }
        }
    s << ")";
    return s;
    }

void
read(std::istream & s, QN & q)
    {
    for(auto& v : q.store()) itensor::read(s,v);
    }

void
write(std::ostream & s, QN const& q)
    {
    for(auto& v : q.store()) itensor::write(s,v);
    }
} //namespace itensor
