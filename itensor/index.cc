//
// Distributed under the ITensor Library License, Version 1.1.
//    (See accompanying LICENSE file.)
//
#include "index.h"

namespace itensor {

using std::string;
using std::stringstream;

std::ostream& 
operator<<(std::ostream& s, const IndexType& it)
    { 
    if(it == Link)      s << "Link"; 
    else if(it == Site) s << "Site"; 
    else if(it == All)  s << "All"; 
    return s; 
    }

int 
IndexTypeToInt(IndexType it)
    {
    if(it == Link) return 1;
    if(it == Site) return 2;
    if(it == All)  return 3;
    Error("No integer value defined for IndexType.");
    return -1;
    }

IndexType 
IntToIndexType(int i)
    {
    if(i == 1) return Link;
    if(i == 2) return Site;
    if(i == 3) return All;
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
        str << "[" << plev << "'s]";
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

string 
nameindex(IndexType it, int plev)
    { 
    static const array<string,3>
    indextypename = {{ "Link","Site", "All" }};
#ifdef DEBUG
    return putprimes(indextypename.at(int(it)),plev);
#else
    return putprimes(indextypename[int(it)],plev); 
#endif
    }

string 
nameint(const string& f, int n)
    { 
    stringstream ss; 
    ss << f << n; 
    return ss.str(); 
    }


//
// IndexDat
// Storage for Index objects.
//
struct IndexDat
    {
    //////////////

    typedef mt19937 IDGenerator;
    typedef IDGenerator::result_type IDType;

    const IDType id;
    const int m;
    const IndexType type;
    const string sname;

    //////////////

    IndexDat(const string& ss, int mm, IndexType it, IDType id);

    static const IndexDatPtr&
    Null();

    private:

    //These methods are not implemented
    //to disallow copying
    IndexDat(const IndexDat&);
    void operator=(const IndexDat&);

    }; //class IndexDat

IndexDat::
IndexDat(const string& ss, int m_, IndexType it, IDType id_)
    : 
    id(id_),
    m(m_), 
    type(it), 
    sname(ss)
    { }

const IndexDatPtr& IndexDat::
Null()
    {
    static IndexDatPtr Null_ = itensor::make_shared<IndexDat>("Null",1,Site,0);
    return Null_;
    }

//
// class Index
//


IndexDat::IDType 
generateID()
    {
    static IndexDat::IDGenerator rng(std::time(NULL) + getpid());
    return rng();

    //static IDType nextid = 0;
    //++nextid;
    //return nextid;
    }


Index::
Index() 
    : 
    p(IndexDat::Null()), 
    primelevel_(0) 
    { 
    //setUniqueReal();
    }

Index::
Index(const string& name, int mm, IndexType it, int plev) 
    : 
    p(itensor::make_shared<IndexDat>(name,mm,it,generateID())), 
    primelevel_(plev) 
    { 
    if(it == All) Error("Constructing Index with type All disallowed");
    //setUniqueReal();
    }


int Index::
m() const { return p->m; }

IndexType Index::
type() const { return p->type; }

string Index::
name() const  { return putprimes(rawname(),primelevel_); }

const string& Index::
rawname() const 
    { 
    return p->sname; 
    }

bool Index::
valid() const { return (p != IndexDat::Null()); }

int Index::
primeLevel() const { return primelevel_; }

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

bool Index::
operator==(const Index& other) const 
    { 
    return (p->id == other.p->id) && (primelevel_ == other.primelevel_); 
    }

bool Index::
noprimeEquals(const Index& other) const
    { 
    return (p->id == other.p->id);
    }

bool Index::
operator<(const Index& other) const 
    { 
    return (p->id < other.p->id);
    }

IndexVal Index::
operator()(int i) const { return IndexVal(*this,i); }

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

    const int t = IndexTypeToInt(p->type);
    s.write((char*) &t,sizeof(t));

    s.write((char*) &(p->id),sizeof(p->id));

    s.write((char*) &(p->m),sizeof(p->m));

    const int nlength = p->sname.length();
    s.write((char*) &nlength,sizeof(nlength));

    s.write(p->sname.data(),nlength+1);
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

    IndexDat::IDType id;
    s.read((char*) &id, sizeof(id));

    int mm; 
    s.read((char*) &mm,sizeof(mm));

    int nlength; 
    s.read((char*) &nlength,sizeof(nlength));

    char* newname = new char[nlength+1]; 
    s.read(newname,nlength+1);
    string ss(newname); 
    delete newname;

    p = itensor::make_shared<IndexDat>(ss,mm,IntToIndexType(t),id);
    return *this;
    }

const Index& Index::
Null()
    {
    static const Index Null_;
    return Null_;
    }

std::ostream& 
operator<<(std::ostream& s, const Index& t)
    {
    if(t.name() != "" && t.name() != " ") s << t.name();
    return s << "(" << nameindex(t.type(),t.primeLevel()) 
             << "," << (t.p->id % 10000) << "):" << t.m();
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

}; //namespace itensor

