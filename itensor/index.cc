//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "index.h"

namespace itensor {

using std::string;
using std::stringstream;

const char*
indexTypeName(IndexType it)
    {
    if(it == Link)      return "Link"; 
    else if(it == Site) return "Site"; 
    else if(it == All)  return "All"; 
    else if(it == NullIndex)  return "NullIndex"; 
    else return "";
    }

std::ostream& 
operator<<(std::ostream& s, IndexType it)
    { 
    s << indexTypeName(it);
    return s; 
    }

int 
IndexTypeToInt(IndexType it)
    {
    if(it == Link) return 1;
    if(it == Site) return 2;
    if(it == All)  return 3;
    if(it == NullIndex) return 4;
    Error("No integer value defined for IndexType.");
    return -1;
    }

IndexType 
IntToIndexType(int i)
    {
    if(i == 1) return Link;
    if(i == 2) return Site;
    if(i == 3) return All;
    if(i == 4) return NullIndex;
    printfln("No IndexType value defined for i=%d\n",i);
    Error("Undefined IntToIndexType value");
    return Link;
    }

string 
putprimes(string s, int plev)
    { 
    stringstream str;
    str << s;
    if(plev < 0) Error("Negative prime level");
    if(plev > 3)
        {
        str << "[" << plev << "']";
        }
    else
        {
        for(int i = 1; i <= plev; ++i) str << "\'";
        }
    return str.str();
    }

string 
nameindex(IndexType it, int plev)
    { 
    return putprimes(indexTypeName(it),plev); 
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
    if(!valid()) Error("Index::write: Index is default initialized");

    s.write((char*) &primelevel_,sizeof(primelevel_));

    const auto t = IndexTypeToInt(type_);
    s.write((char*) &t,sizeof(t));

    s.write((char*) &(id_),sizeof(id_));

    s.write((char*) &(m_),sizeof(m_));

    const int nlength = sname_.length();
    s.write((char*) &nlength,sizeof(nlength));

    s.write(sname_.c_str(),nlength+1);
    }

Index& Index::
read(std::istream& s)
    {
    s.read((char*) &primelevel_,sizeof(primelevel_));
#ifdef DEBUG
    if(primelevel_ < 0)
        {
        Error("Negative primeLevel");
        }
#endif

    int t; 
    s.read((char*) &t,sizeof(t));
    type_ = IntToIndexType(t);

    s.read((char*) &id_, sizeof(id_));

    s.read((char*) &m_,sizeof(m_));

    int nlength; 
    s.read((char*) &nlength,sizeof(nlength));

    auto newname = std::unique_ptr<char[]>(new char[nlength+1]);
    s.read(newname.get(),nlength+1);
    sname_ = string(newname.get()); 

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
    i(0) 
    { }

IndexVal::
IndexVal(const Index& index_, long i_) 
    : 
    index(index_),
    i(i_)
    { 
#ifdef DEBUG
    if(!index) Error("IndexVal initialized with default initialized Index");
    if(i_ < 1 || i_ > index.m())
        {
        println("i = ",i_);
        println("index = ",index);
        Error("i out of range");
        }
#endif
    }

bool IndexVal::
operator==(const IndexVal& other) const
    {
    return (index == other.index && i == other.i);
    }


string
showm(const Index& I) { return nameint("m=",I.m()); }

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv)
    { 
    const Index& ii = iv.index;
    return s << "IndexVal: i = " << iv.i << ", ind = " << ii << "\n"; 
    }

}; //namespace itensor

