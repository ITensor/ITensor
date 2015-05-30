//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <array>
#include "itensor/index.h"
#include "itensor/util/readwrite.h"

namespace itensor {

using std::string;
using std::stringstream;


string 
putprimes(string s, int plev)
    { 
    stringstream str;
    str << s;
    if(plev < 0) Error("Negative prime level");
    if(plev > 3)
        {
        //str << "[" << plev << "']";
        str << "^" << plev;
        }
    else
        {
        for(int i = 1; i <= plev; ++i) str << "\'";
        }
    return str.str();
    }

string 
nameindex(const IndexType& it, int plev)
    { 
    return putprimes(it.c_str(),plev); 
    }

string 
nameint(const string& f, int n)
    { 
    stringstream ss; 
    ss << f << n; 
    return ss.str(); 
    }

//
// class Index
//


string Index::
name() const  { return putprimes(sname_,primelevel_); }

Index& Index::
primeLevel(int plev) 
    { 
    primelevel_ = plev; 
#ifdef DEBUG
    if(primelevel_ < 0)
        Error("Negative primeLevel");
#endif
    return *this;
    }


Index& Index::
mapprime(int plevold, int plevnew, IndexType type)
    {
    if(primelevel_ == plevold)
        {
        if(type == All || type == this->type())
            {
            primelevel_ = plevnew;
#ifdef DEBUG
            if(primelevel_ < 0)
                {
                Error("Negative primeLevel");
                }
#endif
            }
        }
    return *this;
    }


Index& Index::
prime(IndexType type, int inc)
    {
    if(type == this->type() || type == All)
        {
        primelevel_ += inc;
#ifdef DEBUG
        if(primelevel_ < 0)
            {
            Error("Increment led to negative primeLevel");
            }
#endif
        }
    return *this;
    }

void Index::
write(std::ostream& s) const 
    { 
    if(!bool(*this)) Error("Index::write: Index is default initialized");
    itensor::write(s,primelevel_);
    itensor::write(s,type_);
    itensor::write(s,id_);
    itensor::write(s,m_);
    itensor::write(s,sname_);
    }

Index& Index::
read(std::istream& s)
    {
    itensor::read(s,primelevel_);
    itensor::read(s,type_);
    itensor::read(s,id_);
    itensor::read(s,m_);
    itensor::read(s,sname_);

#ifdef DEBUG
    if(primelevel_ < 0) Error("Negative primeLevel");
#endif

    return *this;
    }

std::ostream& 
operator<<(std::ostream& s, const Index& t)
    {
    if(t.name() != "" && t.name() != " ") s << t.name();
    return s << "(" << nameindex(t.type(),t.primeLevel()) 
             << "," << (t.id() % 10000) << "):" << t.m();
    }

IndexVal::
IndexVal() 
    : 
    val(0) 
    { }

IndexVal::
IndexVal(const Index& index_, long val_) 
    : 
    index(index_),
    val(val_)
    { 
#ifdef DEBUG
    if(!index) Error("IndexVal initialized with default initialized Index");
    //Can also use IndexVal's to indicate prime increments:
    //if(val_ < 1 || val_ > index.m())
    //    {
    //    println("val = ",val_);
    //    println("index = ",index);
    //    Error("val out of range");
    //    }
#endif
    }

bool IndexVal::
operator==(const IndexVal& other) const
    {
    return (index == other.index && val == other.val);
    }


string
showm(const Index& I) { return nameint("m=",I.m()); }

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv)
    { 
    const Index& ii = iv.index;
    return s << "IndexVal: val = " << iv.val << ", ind = " << ii << "\n"; 
    }

void
add(Args& args, 
    const Args::Name& name, 
    IndexType it) 
    { 
    args.add(name,it.c_str()); 
    }

IndexType
getIndexType(const Args& args, 
             const Args::Name& name)
    {
    if(!args.defined(name)) Error(format("Name %s not found in Args",name));
    return IndexType(args.getString(name).c_str());
    }

IndexType
getIndexType(const Args& args, 
             const Args::Name& name, 
             IndexType default_val)
    {
    if(!args.defined(name)) return default_val; 
    return IndexType(args.getString(name).c_str());
    }

} //namespace itensor

