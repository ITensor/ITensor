//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
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
        str << "'" << plev;
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
nameint(string const& f, int n)
    { 
    return format("%s%d",f,n);
    }

//
// class Index
//



Index::id_type Index::
generateID()
    {
    static thread_local Index::IDGenerator G;
    return G();
    }

Index::
Index() 
    : 
    id_(0),
    primelevel_(0),
    m_(1),
    type_(NullInd)
    { }

Index::
Index(const std::string& name, long m, IndexType type, int plev) 
    : 
    id_(generateID()),
    primelevel_(plev),
    m_(m),
    type_(type),
    name_(name.c_str())
    { 
#ifdef DEBUG
    if(type_ == All) Error("Constructing Index with type All disallowed");
    if(type_ == NullInd) Error("Constructing Index with type NullInd disallowed");
#endif
    }


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


string Index::
name() const  { return putprimes(name_.c_str(),primelevel_); }

Index::
operator bool() const { return (id_!=0); }


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


bool Index::
noprimeEquals(Index const& other) const
    { 
    return (id_ == other.id_);
    }

IndexVal Index::
operator()(long val) const
    {
    return IndexVal(*this,val);
    }

Index Index::
operator[](int plev) const
    { 
    auto I = *this;
    I.primeLevel(plev); 
    return I; 
    }

void Index::
write(std::ostream& s) const 
    { 
    if(!bool(*this)) Error("Index::write: Index is default initialized");
    itensor::write(s,primelevel_);
    itensor::write(s,type_);
    itensor::write(s,id_);
    itensor::write(s,m_);
    itensor::write(s,name_);
    }

Index& Index::
read(std::istream& s)
    {
    itensor::read(s,primelevel_);
    itensor::read(s,type_);
    if(Global::read32BitIDs())
        {
        using ID32 = std::mt19937::result_type;
        ID32 oldid = 0;
        itensor::read(s,oldid);
        id_ = oldid;
        }
    else
        {
        itensor::read(s,id_);
        }
    itensor::read(s,m_);
    itensor::read(s,name_);

#ifdef DEBUG
    if(primelevel_ < 0) Error("Negative primeLevel");
#endif

    return *this;
    }

bool 
operator==(Index const& i1, Index const& i2)
    { 
    return (i1.id() == i2.id()) && (i1.primeLevel() == i2.primeLevel()); 
    }

bool 
operator!=(Index const& i1, Index const& i2)
    { 
    return not operator==(i1,i2);
    }

bool
operator>(Index const& i1, Index const& i2)
    { 
    if(i1.m() == i2.m()) 
        {
        if(i1.id() == i2.id()) return i1.primeLevel() > i2.primeLevel();
        return i1.id() > i2.id();
        }
    return i1.m() > i2.m();
    }

bool
operator<(Index const& i1, Index const& i2)
    {
    if(i1.m() == i2.m()) 
        {
        if(i1.id() == i2.id()) return i1.primeLevel() < i2.primeLevel();
        return i1.id() < i2.id();
        }
    return i1.m() < i2.m();
    }




std::ostream& 
operator<<(std::ostream & s, Index const& t)
    {
    s << "(\"" << t.rawname();
    s << "\"," << t.m();
    s << "," << t.type().c_str();
    if(Global::showIDs()) 
        {
        s << "|" << (t.id() % 1000);
        //s << "," << t.id();
        }
    s << ")"; 
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
    return s;
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

bool
operator==(IndexVal const& iv1, IndexVal const& iv2)
    {
    return (iv1.index == iv2.index && iv1.val == iv2.val);
    }

bool
operator!=(IndexVal const& iv1, IndexVal const& iv2)
    {
    return not operator==(iv1,iv2);
    }

bool
operator==(Index const& I, IndexVal const& iv)
    {
    return iv.index == I;
    }

bool
operator==(IndexVal const& iv, Index const& I)
    {
    return iv.index == I;
    }

Index
sim(Index const& I, int plev)
    {
    return Index("~"+I.rawname(),I.m(),I.type(),plev);
    }

string
showm(Index const& I) { return nameint("m=",I.m()); }


std::ostream& 
operator<<(std::ostream& s, IndexVal const& iv)
    { 
    return s << "IndexVal: val = " << iv.val 
             << ", ind = " << iv.index;
    }

void
add(Args            & args, 
    Args::Name const& name, 
    IndexType         it) 
    { 
    args.add(name,it.c_str()); 
    }

IndexType
getIndexType(Args       const& args, 
             Args::Name const& name)
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

