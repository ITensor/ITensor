#include "itensor/qn.h"
#include "itensor/util/readwrite.h"
#include "itensor/util/print_macro.h"

namespace itensor {

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


QNVal&
operator+=(QNVal& qva, QNVal const& qvb) 
    { 
    checkCompatible(qva,qvb);
    qva.set(qva.val()+qvb.val());
    return qva;
    }

QNVal&
operator-=(QNVal& qva, QNVal const& qvb) 
    { 
    checkCompatible(qva,qvb);
    qva.set(qva.val()-qvb.val());
    return qva;
    }

QNVal&
operator*=(QNVal& qva, Arrow dir)
    { 
    qva.set(qva.val() * static_cast<int>(dir));
    return qva;
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

std::ostream& 
operator<<(std::ostream & s, QNVal const& qv)
    {
    s << "{";
    if(qv.name() != QNName())
        {
        s << "\"" << qv.name() << "\",";
        }
    s << qv.val() << "," << qv.mod() << "}";
    return s;
    }

////////////////////////
////////////////////////
//
// QN
//
////////////////////////
////////////////////////

QN::
QN(qn_t q0)
    {
    qvs_[0] = QNVal(q0);
    }

QN::
QN(qn_t q0,
   qn_t q1)
    {
    addVal(QNVal(q0));
    addVal(QNVal(q1));
    }

QN::
QN(qn_t q0,
   qn_t q1,
   qn_t q2)
    {
    addVal(QNVal(q0));
    addVal(QNVal(q1));
    addVal(QNVal(q2));
    }

QN::
QN(qn_t q0,
   qn_t q1,
   qn_t q2,
   qn_t q3)
    {
    addVal(QNVal(q0));
    addVal(QNVal(q1));
    addVal(QNVal(q2));
    addVal(QNVal(q3));
    }

//template<typename Func>
//void
//mergeDo(QN & q, 
//        QNVal const& v,
//        Func && eq_action)
//    {
//    auto& qvs = q.store();
//    auto n = QNSize()-1;
//    while(n >= 1)
//        {
//        Print(n);
//        if(isActive(qvs[n-1]))
//            {
//            if(qvs[n-1].name() == v.name()) 
//                {
//                eq_action(qvs[n-1],v);
//                return;
//                }
//            else if(qvs[n-1].name() < v.name())
//                {
//                qvs[n] = v;
//                return;
//                }
//            else //n-1's name is > v.name()
//                {
//                //move qvs_[n-1] to the right
//                //and leave a 'hole' at n
//                qvs[n] = qvs[n-1]; 
//                }
//            }
//        n -= 1;
//        }
//    assert(n == 0);
//    qvs[n] = v;
//    }


void QN::
addVal(QNVal const& qv)
    {
    if(isActive(qvs_.back())) Error("addVal: all QN slots are filled");

    auto n = QNSize()-1;
    while(n >= 1)
        {
        if(isActive(qvs_[n-1]))
            {
            if(qvs_[n-1].name() == qv.name()) 
                {
                Print(qvs_[n-1].name());
                Print(qv.name());
                Error("Duplicate name in QN");
                }
            else if(qvs_[n-1].name() < qv.name())
                {
                qvs_[n] = qv;
                return;
                }
            else //n-1's name is > qv.name()
                {
                //move qvs_[n-1] to the right
                //and leave a 'hole' at n
                qvs_[n] = qvs_[n-1]; 
                }
            }
        n -= 1;
        }
    assert(n == 0);
    qvs_[n] = qv;
    }

QN::qn_t QN::
getVal(QNName const& name) const
    {
    for(auto& v : qvs_)
        {
        if(v.name() == name) return v.val();
        }
    return 0;
    }

void QN::
modAssign(QN const& qo)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        qvs_[n] = QNVal(qvs_[n].name(),qvs_[n].val(),qo.qvs_[n].mod());
        }
    }

void
printFull(QN const& q)
    {
    print("QN(");
    for(auto n : range(q.store()))
        {
        auto& v = q.store()[n];
        if(!isActive(v)) break;
        if(n > 0) print(",");
        print("{\"",v.name(),"\",",v.val(),",",v.mod(),"}");
        }
    println(")");
    }

void
checkCompatible(QN const& qa, QN const& qb) 
    {
#ifdef DEBUG
    //for(auto n : range(QNSize()))
    //    {
    //    }
#endif
    }

bool
operator==(QN qa, QN const& qb)
    {
    for(auto& bv : qb.store()) if(bv.val() != 0)
        {
        for(auto& av : qa.store())
            {
            if(av.name() == bv.name() 
               && av.val() != bv.val()) return false;
            }
        }
    for(auto& av : qa.store()) if(av.val() != 0)
        {
        for(auto& bv : qb.store())
            {
            if(bv.name() == av.name() 
               && bv.val() != av.val()) return false;
            }
        }
    return true;
    }

bool
operator<(QN const& qa, QN const& qb)
    {
    for(auto n : range(QNSize()))
        {
        auto& vn = qa.val0(n);
        auto bval = qb.getVal(vn.name());
        if(vn.val() == bval) continue;
        return vn.val() < bval;
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

//
// Common implementation for += and -= of QNs
//
template<typename Func>
void
combineQN(QN & qa, 
          QN const& qb,
          Func && operation)
    { 
    auto& avs = qa.store();
    for(auto & bv : qb.store())
        {
        if(not isActive(bv)) break;

        for(auto n : range(QNSize()))
            {
            auto& av = avs[n];
            if(av.name() == bv.name()) 
                {
                operation(av,bv);
                break;
                }
            else if(not isActive(avs[n]))
                {
                av = bv;
                break;
                }
            else if(bv.name() < av.name() && 
                    (n==0 || bv.name() > avs[n-1].name()))
                {
                //Name of bv is missing from qa,
                //move this and all remaining vals over one place
                for(size_t m = QNSize()-1; m > n; m -= 1)
                    {
                    avs[m] = avs[m-1];
                    }
                //Put bv into n'th place
                av = bv;
                break;
                }
            }
        }
    }

QN&
operator+=(QN & qa, QN const& qb) 
    { 
    auto addOp = [](QNVal & qva, QNVal const& qvb)
        {
        qva += qvb;
        };
    combineQN(qa,qb,addOp);
    return qa;
    }

QN&
operator-=(QN & qa, QN const& qb) 
    { 
    qa += (-qb);
    return qa;
    }

QN&
operator*=(QN & qa, Arrow dir)
    { 
    for(auto& v : qa.store())
        {
        v *= dir;
        }
    return qa;
    }

std::ostream& 
operator<<(std::ostream & s, QN const& q)
    {
    s << "QN(";
    for(auto n : range(q.store()))
        {
        auto& v = q.store()[n];
        if(!isActive(v)) break;
        if(n > 0) s << ",";
        if(v.mod() == 1)
            {
            s << "{\"" << v.name() << "\"," << v.val() << "}";
            }
        else
            {
            s << "{\"" << v.name() << "\"," << v.val() << "," << v.mod() << "}";
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
