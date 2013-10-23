//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "index.h"
#include "boost/make_shared.hpp"
#include "boost/random/lagged_fibonacci.hpp"

using namespace std;
using boost::format;
//using boost::shared_ptr;
//using boost::make_shared;


ostream& 
operator<<(ostream& s, const IndexType& it)
    { 
    if(it == Link) s << "Link"; 
    else if(it == Site) s << "Site"; 
    else if(it == All) s << "All"; 
    return s; 
    }

int 
IndexTypeToInt(IndexType it)
    {
    if(it == Link) return 1;
    if(it == Site) return 2;
    if(it == All) return 3;
    Error("No integer value defined for IndexType.");
    return -1;
    }

IndexType 
IntToIndexType(int i)
    {
    if(i == 1) return Link;
    if(i == 2) return Site;
    if(i == 3) return All;
    cout << format("No IndexType value defined for i=%d\n")%i 
              << endl;
    Error("Undefined IntToIndexType value");
    return Link;
    }

string 
putprimes(string s, int plev)
    { 
    for(int i = 1; i <= plev; ++i) 
        s += "\'"; 
    return s;
    }

string 
nameindex(IndexType it, int plev)
    { 
    static const Array<string,3>
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
    //
    // Public Data Members

    const IndexType type;
    const int m;
    const Real ur;
    const string sname;

    //
    //////////////

    IndexDat(const string& ss, int mm, IndexType it, Real ur);

    static const IndexDatPtr&
    Null();

    private:

    //These methods are not implemented
    //to disallow copying
    IndexDat(const IndexDat&);
    void operator=(const IndexDat&);

    }; //class IndexDat

IndexDat::
IndexDat(const string& ss, int m_, IndexType it, Real ur_)
    : 
    type(it), 
    m(m_), 
    ur(ur_),
    sname(ss)
    { }

const IndexDatPtr& IndexDat::
Null()
    {
    static IndexDatPtr Null_ = boost::make_shared<IndexDat>("Null",1,Site,0.0);
    return Null_;
    }

//
// class Index
//

//typedef boost::random::lagged_fibonacci1279 
typedef boost::random::lagged_fibonacci2281 
Generator;

Real 
generateUniqueReal()
    {
    static const char seed = 's';

    //Construct rng and seed with address of seed
    static Generator rng((uintptr_t)&seed);

    return rng();
    }



Index::
Index() 
    : 
    p(IndexDat::Null()), 
    primelevel_(0) 
    { }

Index::
Index(const string& name, int mm, IndexType it, int plev) 
    : 
    p(boost::make_shared<IndexDat>(name,mm,it,generateUniqueReal())), 
    primelevel_(plev) 
    { 
    if(it == All) Error("Constructing Index with type All disallowed");
    }

Index::
Index(const IndexDatPtr& p_, int plev) 
    : 
    p(p_),
    primelevel_(plev) 
    { }

int Index::
m() const { return p->m; }

IndexType Index::
type() const { return p->type; }

string Index::
name() const  { return putprimes(rawname(),primelevel_); }

const string& Index::
rawname() const { return p->sname; }


//static const Real real_min = std::numeric_limits<Real>::min();

Real Index::
uniqueReal() const { return p->ur*(1+(primelevel_/10.)); }

bool Index::
isNull() const { return (p == IndexDat::Null()); }

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
    return (p->ur == other.p->ur && primelevel_ == other.primelevel_); 
    }

bool Index::
noprimeEquals(const Index& other) const
    { 
    return (p->ur == other.p->ur); 
    }

bool Index::
operator<(const Index& other) const 
    { return (uniqueReal() < other.uniqueReal()); }

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
write(ostream& s) const 
    { 
    if(isNull()) Error("Index::write: Index is null");

    s.write((char*) &primelevel_,sizeof(primelevel_));

    const int t = IndexTypeToInt(p->type);
    s.write((char*) &t,sizeof(t));

    s.write((char*) &(p->ur),sizeof(p->ur));

    s.write((char*) &(p->m),sizeof(p->m));

    const int nlength = p->sname.length();
    s.write((char*) &nlength,sizeof(nlength));

    s.write(p->sname.data(),nlength+1);
    }

Index& Index::
read(istream& s)
    {
    s.read((char*) &primelevel_,sizeof(primelevel_));
#ifdef DEBUG
    if(primelevel_ < 0)
        {
        Error("Negative primeLevel");
        }
#endif

    int t; s.read((char*) &t,sizeof(t));

    Real ur;
    s.read((char*) &ur, sizeof(ur));

    int mm; 
    s.read((char*) &mm,sizeof(mm));

    int nlength; 
    s.read((char*) &nlength,sizeof(nlength));

    char* newname = new char[nlength+1]; 
    s.read(newname,nlength+1);
    string ss(newname); 
    delete newname;

    p = boost::make_shared<IndexDat>(ss,mm,IntToIndexType(t),ur);

    return *this;
    }

const Index& Index::
Null()
    {
    static const Index Null_;
    return Null_;
    }

ostream& 
operator<<(ostream& s, const Index& t)
    {
    if(t.name() != "" && t.name() != " ") s << t.name();
    const int iur = (int) abs(10000*deprimed(t).uniqueReal());
    return s << "(" << nameindex(t.type(),t.primeLevel()) 
             << "," << iur << "):" << t.m();
    }

IndexVal::
IndexVal() 
    : 
    i(0) 
    { }

IndexVal::
IndexVal(const Index& index, int i_) 
    : 
    Index(index),
    i(i_)
    { 
#ifdef DEBUG
    if(index == Index::Null())
        {
        Error("IndexVal initialized with null Index");
        }
    if(i_ < 1 || i_ > index.m())
        {
        cout << "i = " << i_ << endl;
        cout << "index = " << index << endl;
        Error("i out of range");
        }
#endif
    }


const IndexVal& IndexVal::
Null()
    {
    static const IndexVal Null_;
    return Null_;
    }

string
showm(const Index& I) { return nameint("m=",I.m()); }

ostream& 
operator<<(ostream& s, const IndexVal& iv)
    { 
    const Index& ii = iv;
    return s << "IndexVal: i = " << iv.i << ", ind = " << ii << "\n"; 
    }

