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

//string 
//nameindex(const IndexType& it, int plev)
//    { 
//    return putprimes(it.c_str(),plev); 
//    }

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
    m_(1),
    primelevel_(0),
    tags_(TagSet())
    {
    }

Index::
Index(long m, const TagSet& t, int plev)
    :
    id_(generateID()),
    m_(m),
    primelevel_(plev),
    tags_(t)
    { 
    } 

Index::
Index(long m, int plev)
    :
    id_(generateID()),
    m_(m),
    primelevel_(plev),
    tags_(TagSet())
    { 
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


Index::
operator bool() const { return (id_!=0); }


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

//TODO: read and write for TagSet
void Index::
write(std::ostream& s) const 
    { 
    if(!bool(*this)) Error("Index::write: Index is default initialized");
    itensor::write(s,primelevel_);
    //itensor::write(s,tags_);
    itensor::write(s,id_);
    itensor::write(s,m_);
    }

//TODO: read and write for TagSet
Index& Index::
read(std::istream& s)
    {
    itensor::read(s,primelevel_);
    //itensor::read(s,tags_);
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

#ifdef DEBUG
    if(primelevel_ < 0) Error("Negative primeLevel");
#endif

    return *this;
    }

bool 
operator==(Index const& i1, Index const& i2)
    { 
    return (i1.id() == i2.id()) && (tags(i1) == tags(i2)) && (i1.primeLevel() == i2.primeLevel()); 
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
    s << "(" << t.m();
    if(size(tags(t)) > 0) s << "," << tags(t);
    //s << "," << t.type().c_str();
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
    return Index(I.m(),I.tags(),plev);
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
    TagSet            ts) 
    { 
    args.add(name,ts.c_str()); 
    }

TagSet
getTagSet(Args       const& args, 
          Args::Name const& name)
    {
    if(!args.defined(name)) Error(format("Name %s not found in Args",name));
    return TagSet(args.getString(name).c_str());
    }

TagSet
getTagSet(const Args& args, 
          const Args::Name& name, 
          TagSet default_val)
    {
    if(!args.defined(name)) return default_val; 
    return TagSet(args.getString(name).c_str());
    }

//void
//add(Args            & args, 
//    Args::Name const& name, 
//    TagSet            it) 
//    { 
//    args.add(name,it.c_str()); 
//    }

//IndexType
//getIndexType(Args       const& args, 
//             Args::Name const& name)
//    {
//    if(!args.defined(name)) Error(format("Name %s not found in Args",name));
//    return IndexType(args.getString(name).c_str());
//    }

//IndexType
//getIndexType(const Args& args, 
//             const Args::Name& name, 
//             IndexType default_val)
//    {
//    if(!args.defined(name)) return default_val; 
//    return IndexType(args.getString(name).c_str());
//    }

} //namespace itensor

