//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "index.h"

using namespace std;
using boost::format;
using boost::shared_ptr;
using boost::make_shared;

struct UniqueID
    {
    boost::uuids::uuid id;

    UniqueID() : id(boost::uuids::random_generator()()) { }

    UniqueID& operator++();

    operator boost::uuids::uuid() const { return id; }

    friend std::ostream&
    operator<<(std::ostream& s, const UniqueID& uid);
    };

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


//
// IndexDat
// Storage for Index objects.
//
class IndexDat
    {
    public:

    //////////////
    //
    // Public Data Members

    const IndexType _type;
    const boost::uuids::uuid ind;
    const int m_;
    const Real ur;
    const std::string sname;

    //
    //////////////

    IndexDat(const std::string& name="", int mm = 1,IndexType it=Link);

    //For use with read/write functionality of Index class
    IndexDat(const std::string& ss, int mm, IndexType it, const boost::uuids::uuid& ind_);

    explicit
    IndexDat(Index::Imaker im);

    static const IndexDatPtr&
    Null();

    static const IndexDatPtr&
    ReImDat();

    static const IndexDatPtr&
    ReImDatP();

    static const IndexDatPtr&
    ReImDatPP();

    private:

    //These constructors are not implemented
    //to disallow copying
    IndexDat(const IndexDat&);
    void operator=(const IndexDat&);

    }; //class IndexDat



const UniqueID& //IndexDat::
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

Real 
setUniqueReal(const boost::uuids::uuid ind, IndexType type)
    {
    //return sin(ind * sqrt(1.0/7.0) + ((int)type - (int)Site) * sqrt(1.0 / 13.0));
    Real arg = 0;
    int pn = 1;
    for(int i = int(ind.size())-1; i >= 0; --i)
        { arg += ind.data[i]*sqrt(1.0/(prime_number(++pn))); }
    arg *= sqrt(1.0/(prime_number(++pn)));
    arg += ((int)type - (int)Site) * sqrt(1.0/(prime_number(++pn)));
    return sin(arg);
    }

IndexDat::
IndexDat(const std::string& name, int mm,IndexType it) 
    : 
    _type(it), 
    ind(nextID()),
    m_(mm), 
    ur(setUniqueReal(ind,it)),
    sname(name)
    { 
    if(it == ReIm) Error("Constructing Index with type ReIm disallowed");
    if(it == All) Error("Constructing Index with type All disallowed");
    }

IndexDat::
IndexDat(const std::string& ss, int mm, IndexType it, const boost::uuids::uuid& ind_)
    : 
    _type(it), 
    ind(ind_), 
    m_(mm), 
    ur(setUniqueReal(ind_,it)),
    sname(ss)
    { 
    if(it == ReIm) Error("Constructing Index with type ReIm disallowed");
    if(it == All) Error("Constructing Index with type All disallowed");
    }

const char*
staticSetName(Index::Imaker im)
    {
    if(im == Index::makeNull)
        return "Null";
    else
    if(im == Index::makeReIm)
        return "ReIm";
    else
    if(im == Index::makeReImP)
        return "ReImP";
    else
    if(im == Index::makeReImPP)
        return "ReImPP";
    Error("Imaker case not handled");
    return "";
    }

boost::uuids::uuid
staticSetInd(Index::Imaker im)
    {
    //Don't use random uuid generator for these static IndexDats
    boost::uuids::string_generator gen;
    if(im==Index::makeNull)
        return gen("{00000000-0000-0000-0000-000000000000}"); 
    else //if im == ReIm, ReImP, or ReImPP                               
        return gen("{10000000-0000-0000-0000-000000000000}"); 
    Error("Imaker case not handled");
    return gen("{00000000-0000-0000-0000-000000000000}");
    }

IndexDat::
IndexDat(Index::Imaker im) 
    : 
    _type((im==Index::makeNull ? Site : ReIm)), 
    ind(staticSetInd(im)),
    m_( (im==Index::makeNull) ? 1 : 2),
    ur(im == Index::makeNull ? 0 : setUniqueReal(ind,_type)),
    sname(staticSetName(im))
    { }

const IndexDatPtr& IndexDat::
Null()
    {
    static IndexDatPtr Null_ = make_shared<IndexDat>(Index::makeNull);
    return Null_;
    }

const IndexDatPtr& IndexDat::
ReImDat()
    {
    static IndexDatPtr ReImDat_ = make_shared<IndexDat>(Index::makeReIm);
    return ReImDat_;
    }

const IndexDatPtr& IndexDat::
ReImDatP()
    {
    static IndexDatPtr ReImDatP_ = make_shared<IndexDat>(Index::makeReImP);
    return ReImDatP_;
    }

const IndexDatPtr& IndexDat::
ReImDatPP()
    {
    static IndexDatPtr ReImDatPP_ = make_shared<IndexDat>(Index::makeReImPP);
    return ReImDatPP_;
    }


Index::
Index() 
    : 
    p(IndexDat::Null()), 
    primelevel_(0) 
    { }

Index::
Index(const std::string& name, int mm, IndexType it, int plev) 
    : 
    p(make_shared<IndexDat>(name,mm,it)), 
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
    prime(type,primeinc);
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

std::string
showm(const Index& I) { return nameint("m=",I.m()); }

Real Index::
uniqueReal() const { return p->ur*(1+0.00398406*primelevel_); }

bool Index::
isNull() const { return (p == IndexDat::Null()); }

bool Index::
isNotNull() const { return (p != IndexDat::Null()); }

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
operator()(int i) const { return IndexVal(*this,i); }

void Index::
mapprime(int plevold, int plevnew, IndexType type)
    {
    if(primelevel_ == plevold)
        {
        if((type == All && this->type() != ReIm) || type == this->type())
            primelevel_ = plevnew;
        }
    }

void Index::
prime(int inc) { primelevel_ += inc; }

void Index::
prime(IndexType type, int inc)
    {
    if(type == this->type() ||
       (type == All && this->type() != ReIm))
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
        p = make_shared<IndexDat>(ss,mm,IntToIndexType(t),ind);
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

std::ostream& 
operator<<(std::ostream& s, const IndexVal& iv)
    { 
    return s << "IndexVal: i = " << iv.i << ", ind = " << iv.ind << "\n"; 
    }

