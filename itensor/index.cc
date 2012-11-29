//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "index.h"

using namespace std;
using boost::format;

ostream& 
operator<<(ostream& s, const IndexType& it)
    { 
    if(it == Link) s << "Link"; 
    else if(it == Site) s << "Site"; 
    else if(it == ReIm) s << "ReIm"; 
    else if(it == All) s << "All"; 
    return s; 
    }

int 
IndexTypeToInt(IndexType it)
    {
    if(it == Link) return 1;
    if(it == Site) return 2;
    if(it == ReIm) return 3;
    if(it == All) return 4;
    Error("No integer value defined for IndexType.");
    return -1;
    }

IndexType 
IntToIndexType(int i)
    {
    if(i == 1) return Link;
    if(i == 2) return Site;
    if(i == 3) return ReIm;
    if(i == 4) return All;
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
    static const boost::array<string,4>
    indextypename = {{ "Link","Site","ReIm", "All" }};
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

#define UID_NUM_PRINT 2
std::ostream& 
operator<<(std::ostream& s, const boost::uuids::uuid& id)
    { 
    s.width(2);
    for(boost::uuids::uuid::size_type i = id.size()-UID_NUM_PRINT; i < id.size(); ++i) 
        {
        s << static_cast<unsigned int>(id.data[i]);
        }
    s.width(0);
    return s; 
    }

UniqueID& UniqueID::
operator++()
    {
    int i = id.size(); 
    while(--i >= 0)
        { 
        if(++id.data[i] == 0) continue; 
        break;
        }
    return *this;
    }

std::ostream&
operator<<(std::ostream& s, const UniqueID& uid) 
    { 
    s << uid.id; 
    return s; 
    }

int 
prime_number(int n)
    {
    static const boost::array<int,54> plist = { { 
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 
    79, 83, 89, 97, 101, 103, 107, 109, 113, 
    127, 131, 137, 139, 149, 151, 157, 163, 
    167, 173, 179, 181, 191, 193, 197, 199, 
    211, 223, 227, 229, 233, 239, 241, 251 
    } };
    return plist.at(n);
    }

void 
intrusive_ptr_add_ref(IndexDat* p) 
    { 
    ++(p->numref); 
    }

void 
intrusive_ptr_release(IndexDat* p) 
    { 
    if(!p->is_static_ && --(p->numref) == 0)
        { 
        delete p; 
        } 
    }

const UniqueID& IndexDat::
nextID()
    {
    static UniqueID lastID_;
    static int count_ = 0;
    //After making so many ID's sequentially,
    //call the random number generator again
    if(++count_ > 1000)
        {
        count_ = 0;
        lastID_ = UniqueID();
        }
    return ++lastID_;
    }

void IndexDat::
setUniqueReal()
    {
    //ur = sin(ind * sqrt(1.0/7.0) + ((int)_type - (int)Site) * sqrt(1.0 / 13.0));
    Real arg = 0;
    int pn = 1;
    for(int i = int(ind.size())-1; i >= 0; --i)
        { arg += ind.data[i]*sqrt(1.0/(prime_number(++pn))); }
    arg *= sqrt(1.0/(prime_number(++pn)));
    arg += ((int)_type - (int)Site) * sqrt(1.0/(prime_number(++pn)));
    ur = sin(arg);
    }

IndexDat::
IndexDat(const std::string& name, int mm,IndexType it) 
    : _type(it), 
      ind(nextID()),
      m_(mm), 
      sname(name),
      numref(0),
      is_static_(false)
    { 
    if(it == ReIm) Error("Constructing Index with type ReIm disallowed");
    if(it == All) Error("Constructing Index with type All disallowed");
    setUniqueReal();
    }

IndexDat::
IndexDat(const std::string& ss, int mm, IndexType it, const boost::uuids::uuid& ind_)
    : _type(it), 
      ind(ind_), 
      m_(mm), 
      sname(ss),
      numref(0), 
      is_static_(false)
    { 
    if(it == ReIm) Error("Constructing Index with type ReIm disallowed");
    if(it == All) Error("Constructing Index with type All disallowed");
    setUniqueReal();
    }

IndexDat::
IndexDat(Index::Imaker im) 
    : 
    _type(ReIm), 
    m_( (im==Index::makeNull) ? 1 : 2),
    numref(0), 
    is_static_(true)
    { 
    //Don't use random uuid generator for these static IndexDats
    boost::uuids::string_generator gen;
    if(im==Index::makeNull)
        { ind = gen("{00000000-0000-0000-0000-000000000000}"); }
    else                               
        { ind = gen("{10000000-0000-0000-0000-000000000000}"); }

    if(im == Index::makeNull)
        {
        sname = "Null";
        _type = Site;
        ur = 0.0;
        return;
        }
    else if(im == Index::makeReIm) sname = "ReIm";
    else if(im == Index::makeReImP) sname = "ReImP";
    else if(im == Index::makeReImPP) sname = "ReImPP";
    setUniqueReal(); 
    }

IndexDat* IndexDat::
Null()
    {
    static IndexDat Null_(Index::makeNull);
    return &Null_;
    }

IndexDat* IndexDat::
ReImDat()
    {
    static IndexDat ReImDat_(Index::makeReIm);
    return &ReImDat_;
    }

IndexDat* IndexDat::
ReImDatP()
    {
    static IndexDat ReImDatP_(Index::makeReImP);
    return &ReImDatP_;
    }

IndexDat* IndexDat::
ReImDatPP()
    {
    static IndexDat ReImDatPP_(Index::makeReImPP);
    return &ReImDatPP_;
    }


Index::
Index() 
    : p(IndexDat::Null()), 
      primelevel_(0) 
    { }

Index::
Index(const std::string& name, int mm, IndexType it, int plev) 
    : 
    p(new IndexDat(name,mm,it)), 
    primelevel_(plev) 
    { 
    }

Index::
Index(Imaker im)
    {
    if(im == makeNull)
        p = IndexDat::Null(), primelevel_ = 0;
    else if(im == makeReIm)
        p = IndexDat::ReImDat(), primelevel_ = 0;
    else if(im == makeReImP)
        p = IndexDat::ReImDatP(),  primelevel_ = 1;
    else if(im == makeReImPP)
        p = IndexDat::ReImDatPP(),  primelevel_ = 2;
    else Error("Unrecognized Imaker type.");
    }

Index::
Index(IndexType type,const Index& other, int primeinc) 
    : p(other.p), 
      primelevel_(other.primelevel_)
    {
    primelevel_ = other.primelevel_;
    doprime(type,primeinc);
    }


int Index::
m() const { return p->m_; }

const boost::uuids::uuid& Index::
Ind() const { return p->ind; }

IndexType Index::
type() const { return p->_type; }

std::string Index::
name() const  { return putprimes(rawname(),primelevel_); }

const std::string& Index::
rawname() const { return p->sname; }

void Index::
setname(const std::string& newname) { p->sname = newname; }

std::string Index::
showm() const { return (boost::format("m=%d")%(p->m_)).str(); }

Real Index::
uniqueReal() const { return p->ur*(1+0.00398406*primelevel_); }

bool Index::
isNull() const { return (p == IndexDat::Null()); }

bool Index::
isNotNull() const { return (p != IndexDat::Null()); }

int Index::
count() const { return p->count(); }

int Index::
primeLevel() const { return primelevel_; }

void Index::
primeLevel(int plev) { primelevel_ = plev; }

bool Index::
operator==(const Index& other) const 
    { 
    return (p->ur == other.p->ur && primelevel_ == other.primelevel_); 
    }

bool Index::
noprime_equals(const Index& other) const
    { 
    return (p->ur == other.p->ur); 
    }

bool Index::
operator<(const Index& other) const 
    { return (uniqueReal() < other.uniqueReal()); }

IndexVal Index::
operator()(int i) const 
    { return IndexVal(*this,i); }

void Index::
mapprime(int plevold, int plevnew, IndexType type)
    {
    if(primelevel_ != plevold) return;
    if(this->type() == ReIm) return;
    else if(type == All || type == this->type() )
        {
        primelevel_ = plevnew;
        }
    }

void Index::
doprime(IndexType type, int inc)
    {
    if(this->type() == ReIm) return;
    if(type == All || type == this->type())
        {
        primelevel_ += inc;
        }
    }

void Index::
write(std::ostream& s) const 
    { 
    if(isNull()) Error("Index::write: Index is null");

    s.write((char*) &primelevel_,sizeof(primelevel_));

    const int t = IndexTypeToInt(p->_type);
    s.write((char*) &t,sizeof(t));

    for(int i = 0; i < int(p->ind.size()); ++i) 
        { const char c = p->ind.data[i] - '0'; s.write(&c,sizeof(c)); }

    s.write((char*) &(p->m_),sizeof(p->m_));

    const int nlength = p->sname.length();
    s.write((char*) &nlength,sizeof(nlength));

    s.write(p->sname.data(),nlength+1);
    }

void Index::
read(std::istream& s)
    {
    s.read((char*) &primelevel_,sizeof(primelevel_));

    int t; s.read((char*) &t,sizeof(t));

    boost::uuids::uuid ind;
    for(int i = 0; i < int(ind.size()); ++i) 
        { char c; s.read(&c,sizeof(c)); ind.data[i] = '0'+c; }

    int mm; s.read((char*) &mm,sizeof(mm));

    int nlength; s.read((char*) &nlength,sizeof(nlength));

    char* newname = new char[nlength+1]; 
    s.read(newname,nlength+1);
    std::string ss(newname); 
    delete newname;

    if(IntToIndexType(t) == ReIm)
        {
        if(primelevel_ == 0) 
            p = IndexDat::ReImDat();
        else if(primelevel_ == 1) 
            p = IndexDat::ReImDatP();
        else if(primelevel_ == 2) 
            p = IndexDat::ReImDatPP();
        else
            Error("Illegal primelevel for Index of ReIm type");
        }
    else
        {
        p = new IndexDat(ss,mm,IntToIndexType(t),ind);
        }
    }

const Index& Index::
Null()
    {
    static const Index Null_(makeNull);
    return Null_;
    }

const Index& Index::
IndReIm()
    {
    static const Index IndReIm_(makeReIm);
    return IndReIm_;
    }

const Index& Index::
IndReImP()
    {
    static const Index IndReImP_(makeReImP);
    return IndReImP_;
    }

const Index& Index::
IndReImPP()
    {
    static const Index IndReImPP_(makeReImPP);
    return IndReImPP_;
    }

std::ostream& 
operator<<(std::ostream & s, const Index & t)
    {
    if(t.name() != "" && t.name() != " ") s << t.name() << "/";
    return s << nameindex(t.type(),t.primelevel_) << "-" << t.Ind() << ":" << t.m();
    }

IndexVal::
IndexVal() 
    : ind(Index::Null()),
      i(0) 
    { }

IndexVal::
IndexVal(const Index& index, int i_) 
    : ind(index),
      i(i_)
    { 
#ifdef DEBUG
    if(index == Index::Null())
        Error("IndexVal initialized with null Index");
#endif
    assert(i <= ind.m()); 
    }

IndexVal::
IndexVal(Index::Imaker im)
    {
    if(im == Index::makeNull)
        {
        ind = Index::Null();
        i = 1;
        }
    else
        {
        Error("Imaker type not supported");
        }
    }

bool IndexVal::
operator==(const IndexVal& other) const 
    { 
    return (ind == other.ind && i == other.i); 
    }

IndexVal
primed(const IndexVal& iv, int inc)
    {
    return IndexVal(primed(iv.ind,inc),iv.i);
    }

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv)
    { 
    return s << "IndexVal: i = " << iv.i << ", ind = " << iv.ind << "\n"; 
    }

