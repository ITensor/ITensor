#include "itensor/qn.h"
#include "itensor/util/readwrite.h"
#include "itensor/util/print_macro.h"

using std::string;
using std::vector;

namespace itensor {

void QNum::
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
checkCompatible(QNum const& qva, QNum const& qvb) 
    {
#ifdef DEBUG
    if(qva.name() != qvb.name())
        {
        printfln("qva.name()=%s, qvb.name()=%s",qva.name(),qvb.name());
        Error("Mismatched QNum names");
        }
    if(qva.mod() != qvb.mod())
        {
        printfln("qva.mod()=%d, qvb.mod()=%d",qva.mod(),qvb.mod());
        Error("Mismatched mod factors");
        }
#endif
    }


QNum&
operator+=(QNum& qva, QNum const& qvb) 
    { 
    checkCompatible(qva,qvb);
    qva.set(qva.val()+qvb.val());
    return qva;
    }

QNum&
operator-=(QNum& qva, QNum const& qvb) 
    { 
    checkCompatible(qva,qvb);
    qva.set(qva.val()-qvb.val());
    return qva;
    }

QNum&
operator*=(QNum& qva, Arrow dir)
    { 
    qva.set(qva.val() * static_cast<int>(dir));
    return qva;
    }

bool
operator==(QNum const& qva, QNum const& qvb)
    {
    checkCompatible(qva,qvb);
    return qva.val() == qvb.val();
    }

bool
operator!=(QNum const& qva, QNum const& qvb)
    {
    checkCompatible(qva,qvb);
    return qva.val() != qvb.val();
    }

void
read(std::istream & s, QNum & q)
    {
    QNName n;
    QNum::qn_t v = 0,
                m = 0;
    itensor::read(s,n);
    itensor::read(s,v);
    itensor::read(s,m);
    q = QNum(n,v,m);
    }

void
write(std::ostream & s, QNum const& q)
    {
    itensor::write(s,q.name());
    itensor::write(s,q.val());
    itensor::write(s,q.mod());
    }

std::ostream& 
operator<<(std::ostream & s, QNum const& qv)
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
    qvs_[0] = QNum(q0);
    }

//QN::
//QN(qn_t q0,
//   qn_t q1)
//    {
//    addNum(QNum(q0));
//    addNum(QNum(q1));
//    }
//
//QN::
//QN(qn_t q0,
//   qn_t q1,
//   qn_t q2)
//    {
//    addNum(QNum(q0));
//    addNum(QNum(q1));
//    addNum(QNum(q2));
//    }
//
//QN::
//QN(qn_t q0,
//   qn_t q1,
//   qn_t q2,
//   qn_t q3)
//    {
//    addNum(QNum(q0));
//    addNum(QNum(q1));
//    addNum(QNum(q2));
//    addNum(QNum(q3));
//    }

//template<typename Func>
//void
//mergeDo(QN & q, 
//        QNum const& v,
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
addNum(QNum const& qv)
    {
    if(isActive(qvs_.back())) Error("addNum: all QN slots are filled");

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
val(QNName const& name) const
    {
    for(auto& v : qvs_)
        {
        if(v.name() == name) return v.val();
        }
    return 0;
    }

QN::qn_t QN::
mod(QNName const& name) const
    {
    for(auto& v : qvs_)
        {
        if(v.name() == name) return v.mod();
        }
    return 0;
    }

bool QN::
hasName(QNName const& name) const
    {
    for(auto& v : qvs_)
        {
        if(v.name() == name) return true;
        }
    return false;
    }

void QN::
modAssign(QN const& qo)
    {
    for(size_t n = 0; n < QNSize(); ++n)
        {
        qvs_[n] = QNum(qvs_[n].name(),qvs_[n].val(),qo.qvs_[n].mod());
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
        bool found = false;
        for(auto& av : qa.store())
            {
            if(av.name() == bv.name())
                {
                if(av.val() != bv.val()) return false;
                found = true;
                break;
                }
            }
        if(not found) return false;
        }
    for(auto& av : qa.store()) if(av.val() != 0)
        {
        bool found = false;
        for(auto& bv : qb.store())
            {
            if(bv.name() == av.name())
                {
                if(bv.val() != av.val()) return false;
                found = true;
                break;
                }
            }
        if(not found) return false;
        }
    return true;
    }

bool
operator<(QN const& qa, QN const& qb)
    {
    size_t a = 1;
    size_t b = 1;
    while(a <= QNSize() && b <= QNSize()
          && (isActive(qa.num(a)) || isActive(qb.num(b))))
        {
        if(not isActive(qa.num(a)))
            {
            if(0 == qb.val(b))
                {
                b += 1;
                continue;
                }
            return 0 < qb.val(b);
            }
        else if(not isActive(qb.num(b)))
            {
            if(qa.val(a) == 0)
                {
                a += 1;
                continue;
                }
            return qa.val(a) < 0;
            }
        else //both active
            {
            auto aname = qa.name(a);
            auto bname = qb.name(b);
            if(aname < bname) //b doesn't have aname, assume bval=0
                {
                if(qa.val(a) == 0) 
                    {
                    a += 1;
                    continue;
                    }
                return qa.val(a) < 0;
                }
            else if(bname < aname) //a doesn't have bname, assume aval=0
                {
                if(qb.val(b) == 0)
                    {
                    b += 1;
                    continue;
                    }
                return 0 < qb.val(b);
                }
            else // bname == aname
                {
                if(qa.val(a) == qb.val(b))
                    {
                    a += 1;
                    b += 1;
                    continue;
                    }
                return qa.val(a) < qb.val(b);
                }
            }
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
            if(not isActive(av))
                {
                av = bv;
                break;
                }
            else if(av.name() == bv.name()) 
                {
                operation(av,bv);
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
    auto addOp = [](QNum & qva, QNum const& qvb)
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
        s << "{";
        if(v.name() != QNName()) s << "\"" << v.name() << "\",";
        s << v.val();
        if(v.mod() != 1)
            {
            s << "," << v.mod();
            }
        s << "}";
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

#ifdef ITENSOR_USE_HDF5

void
h5_write(h5::group parent, string const& name, QN const& q)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type","QN",true);
    h5_write_attribute(g,"version",long(1));
    auto names = vector<string>(QNSize());
    auto vals = vector<long>(QNSize());
    auto mods = vector<long>(QNSize());
    for(auto n : range1(QNSize()))
        {
        names[n-1] = q.name(n);
        vals[n-1] = q.val(n);
        mods[n-1] = q.mod(n);
        }
    h5_write(g,"names",names);
    h5_write(g,"vals",vals);
    h5_write(g,"mods",mods);
    }

void
h5_read(h5::group parent, string const& name, QN & q)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<string>(g,"type");
    if(type != "QN") Error("Group does not contain ITensor data in HDF5 file");
    auto names = h5_read<vector<string>>(g,"names");
    auto vals = h5_read<vector<long>>(g,"vals");
    auto mods = h5_read<vector<long>>(g,"mods");
    q = QN();
    for(auto n : range1(QNSize()))
        {
        q.addNum(QNum(names[n-1],vals[n-1],mods[n-1]));
        }
    }

#endif //ITENSOR_USE_HDF5

} //namespace itensor
