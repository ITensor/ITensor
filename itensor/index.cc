//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "itensor/index.h"

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
        str << "[" << plev << "']";
        }
    else
        {
        for(int i = 1; i <= plev; ++i) 
            {
            str << "\'";
            }
        }
    return str.str();
    }

//string static
//nameindex(IndexType const& it, int plev)
//    { 
//    return putprimes(it.c_str(),plev); 
//    }

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


Index::IDType 
generateID()
    {
    static Index::IDGenerator rng(std::time(NULL) + getpid());
    return rng();

    //static IDType nextid = 0;
    //++nextid;
    //return nextid;
    }

//const IndexDatPtr& IndexDat::
//Null()
//    {
//    static IndexDatPtr Null_ = itensor::make_shared<IndexDat>("Null",1,Site,0);
//    return Null_;
//    }

Index::
Index() 
    : 
    id_(0),
    primelevel_(0),
    m_(1),
    type_(NullInd),
    sname_("Null")
    { 
    }

Index::
Index(const string& name, int m, IndexType type, int plev) 
    : 
    id_(generateID()),
    primelevel_(plev),
    m_(m),
    type_(type),
    sname_(name)
    { 
    if(type_ == All) Error("Constructing Index with type All disallowed");
    }


string Index::
name() const  { return putprimes(rawname(),primelevel_); }

bool Index::
valid() const { return (id_ != 0); }

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
prime(int inc) 
    { 
    primelevel_ += inc; 
#ifdef DEBUG
    if(primelevel_ < 0)
        {
        Error("Negative primeLevel");
        }
#endif
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

    type_.write(s);

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

    type_.read(s);

    s.read((char*) &id_, sizeof(id_));

    s.read((char*) &m_,sizeof(m_));

    int nlength; 
    s.read((char*) &nlength,sizeof(nlength));

    auto newname = std::unique_ptr<char[]>(new char[nlength+1]);
    s.read(newname.get(),nlength+1);
    sname_ = string(newname.get()); 

    return *this;
    }

const Index& Index::
Null()
    {
    static const Index Null_;
    return Null_;
    }

std::string Index::
id() const
    {
    return format("%d",id_);
    }

std::ostream& 
operator<<(std::ostream& s, const Index& t)
    {
    s << "(" << t.rawname() << "," << t.m() << ","
      << t.type().c_str() << ")";
    if(t.primeLevel() > 0) 
      {
      if(t.primeLevel() > 3)
        {
        s << "'" << t.primeLevel();
        }
      else
        {
        for(int n = 1; n <= t.primeLevel(); ++n)
          s << "'";
        }
      }
    //{" << (t.id_ % 1000) << "}";
    return s;

    //s << putprimes("",t.primeLevel()) << "(";
    //s << "(";
    //if(t.name() != "" && t.name() != " ") s << t.name();
    //s << "," << t.m() << "," << t.type().c_str() << ")";
    //s << "{" << (t.id_% 1000) << "}";
    //return s;

    //if(t.name() != "" && t.name() != " ") s << "Index(" << t.name();
    //s << "," << nameindex(t.type(),t.primeLevel()) 
    //      << "," << (t.id_ % 10000) << "):" << t.m();
    //return s;
    }

IndexVal::
IndexVal() 
    : 
    i(0) 
    { }

IndexVal::
IndexVal(const Index& index_, int i_) 
    : 
    index(index_),
    i(i_)
    { 
#ifdef DEBUG
    if(index == Index::Null())
        {
        Error("IndexVal initialized with null Index");
        }
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


const IndexVal& IndexVal::
Null()
    {
    static const IndexVal Null_;
    return Null_;
    }

string
showm(const Index& I) { return nameint("m=",I.m()); }

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv)
    { 
    const Index& ii = iv.index;
    return s << "IndexVal: i = " << iv.i << ", ind = " << ii << "\n"; 
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

