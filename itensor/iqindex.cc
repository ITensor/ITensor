//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqindex.h"

using namespace std;
using boost::format;
using boost::array;
using boost::shared_ptr;
using boost::make_shared;

//
// IQIndexDat
//

class IQIndexDat
    {
    public:

    typedef vector<IndexQN>
    StorageT;

    typedef StorageT::iterator
    iterator;
    typedef StorageT::const_iterator
    const_iterator;


    IQIndexDat() { }

    IQIndexDat(const Index& i1, const QN& q1,
               const Index& i2 = Index::Null(), const QN& q2 = QN(),
               const Index& i3 = Index::Null(), const QN& q3 = QN(),
               const Index& i4 = Index::Null(), const QN& q4 = QN(),
               const Index& i5 = Index::Null(), const QN& q5 = QN(),
               const Index& i6 = Index::Null(), const QN& q6 = QN(),
               const Index& i7 = Index::Null(), const QN& q7 = QN(),
               const Index& i8 = Index::Null(), const QN& q8 = QN());

    explicit
    IQIndexDat(vector<IndexQN>& ind_qn);

    IQIndexDat(istream& s);

    const StorageT&
    indices() const { return iq_; }

    int
    size() { return iq_.size(); }

    const Index&
    index(int i) { return iq_[i-1]; }
    const QN&
    qn(int i) { return iq_[i-1].qn; }

    iterator
    begin() { return iq_.begin(); }
    iterator
    end() { return iq_.end(); }

    const_iterator
    begin() const { return iq_.begin(); }
    const_iterator
    end()   const { return iq_.end(); }

    void 
    write(ostream& s) const;

    void 
    read(istream& s);

    static const IQIndexDatPtr& Null();

    static const IQIndexDatPtr& ReImDat();

    static const IQIndexDatPtr& ReImDatP();

    static const IQIndexDatPtr& ReImDatPP();

    void
    makeCopyOf(const IQIndexDat& other);

    private:

    //////////////////

    StorageT iq_;

    /////////////////

    //Disallow copying using =
    void 
    operator=(const IQIndexDat&);

    };

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3,
           const Index& i4, const QN& q4,
           const Index& i5, const QN& q5,
           const Index& i6, const QN& q6,
           const Index& i7, const QN& q7,
           const Index& i8, const QN& q8)
    {
    iq_.push_back(IndexQN(i1,q1));
    if(i2 != Index::Null())
        iq_.push_back(IndexQN(i2,q2));
    if(i3 != Index::Null())
        iq_.push_back(IndexQN(i3,q3));
    if(i4 != Index::Null())
        iq_.push_back(IndexQN(i4,q4));
    if(i5 != Index::Null())
        iq_.push_back(IndexQN(i5,q5));
    if(i6 != Index::Null())
        iq_.push_back(IndexQN(i6,q6));
    if(i7 != Index::Null())
        iq_.push_back(IndexQN(i7,q7));
    if(i8 != Index::Null())
        iq_.push_back(IndexQN(i8,q8));
    }

IQIndexDat::
IQIndexDat(StorageT& ind_qn)
    { 
    iq_.swap(ind_qn); 
    }

void IQIndexDat::
makeCopyOf(const IQIndexDat& other) 
    { 
    iq_ = other.iq_;
    }

IQIndexDat::
IQIndexDat(istream& s) 
    { read(s); }

void IQIndexDat::
write(ostream& s) const
    {
    size_t size = iq_.size();
    s.write((char*)&size,sizeof(size));
    Foreach(const IndexQN& x, iq_)
        { 
        x.write(s); 
        }
    }

void IQIndexDat::
read(istream& s)
    {
    size_t size; s.read((char*)&size,sizeof(size));
    iq_.resize(size);
    Foreach(IndexQN& x, iq_)
        { 
        x.read(s); 
        }
    }

const IQIndexDatPtr& IQIndexDat::
Null()
    {
    static IQIndexDatPtr Null_ = make_shared<IQIndexDat>(Index::Null(),QN());
    return Null_;
    }

const IQIndexDatPtr& IQIndexDat::
ReImDat()
    {
    static IQIndexDatPtr ReImDat_ = make_shared<IQIndexDat>(Index::IndReIm(),QN());
    return ReImDat_;
    }

const IQIndexDatPtr& IQIndexDat::
ReImDatP()
    {
    static IQIndexDatPtr ReImDatP_ = make_shared<IQIndexDat>(Index::IndReImP(),QN());
    return ReImDatP_;
    }

const IQIndexDatPtr& IQIndexDat::
ReImDatPP()
    {
    static IQIndexDatPtr ReImDatPP_ = make_shared<IQIndexDat>(Index::IndReImPP(),QN());
    return ReImDatPP_;
    }

//
// IQIndex Methods
//

#ifdef DEBUG
#define IQINDEX_CHECK_NULL if(pd == 0) Error("IQIndex is null");
#else
#define IQINDEX_CHECK_NULL
#endif

const IQIndexDat::StorageT& IQIndex::
indices() const 
    { 
    IQINDEX_CHECK_NULL
    return pd->indices();
    }

int IQIndex::
nindex() const 
    { 
    IQINDEX_CHECK_NULL
    return (int) pd->size(); 
    }

const Index& IQIndex::
index(int i) const 
    {
    IQINDEX_CHECK_NULL
#ifdef DEBUG
    if(i > nindex())
        {
        Print(nindex());
        Print(i);
        Error("IQIndex::index arg out of range");
        }
#endif
    return pd->index(i);
    }

const QN& IQIndex::
qn(int i) const 
    {
    IQINDEX_CHECK_NULL
#ifdef DEBUG
    if(i > nindex())
        {
        Print(nindex());
        Print(i);
        Error("IQIndex::qn arg out of range");
        }
#endif
    return pd->qn(i);
    }

IQIndex::
IQIndex() 
    : 
    _dir(Neither)
    { }


IQIndex::
IQIndex(const string& name,
                 IndexType it, 
                 Arrow dir, 
                 int plev) 
    : 
    Index(name,1,it,plev), 
    _dir(dir)
    { }

IQIndex::
IQIndex(const string& name, 
        const Index& i1, const QN& q1, 
        Arrow dir) 
    : 
    Index(name,i1.m(),i1.type(),i1.primeLevel()), 
    _dir(dir), 
    pd(make_shared<IQIndexDat>(i1,q1))
    {
    }

IQIndex::
IQIndex(const string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        Arrow dir) 
    : 
    Index(name,i1.m()+i2.m(),i1.type(),i1.primeLevel()), 
    _dir(dir), 
    pd(make_shared<IQIndexDat>(i1,q1,i2,q2))
    {
    if(i2.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        const Index& i3, const QN& q3,
        Arrow dir) 
    : 
    Index(name,i1.m()+i2.m()+i3.m(),i1.type(),i1.primeLevel()), 
    _dir(dir),
    pd(make_shared<IQIndexDat>(i1,q1,i2,q2,i3,q3))
    {
    if(i2.type() != i1.type() 
    || i3.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        const Index& i3, const QN& q3,
        const Index& i4, const QN& q4,
        Arrow dir) 
    : 
    Index(name,i1.m()+i2.m()+i3.m()+i4.m(),i1.type(),i1.primeLevel()), 
    _dir(dir),
    pd(make_shared<IQIndexDat>(i1,q1,i2,q2,i3,q3,i4,q4))
    {
    if(i2.type() != i1.type() 
    || i3.type() != i1.type()
    || i4.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        const Index& i3, const QN& q3,
        const Index& i4, const QN& q4,
        const Index& i5, const QN& q5,
        Arrow dir) 
    : 
    Index(name,i1.m()+i2.m()+i3.m()+i4.m()+i5.m(),i1.type(),i1.primeLevel()), 
    _dir(dir),
    pd(new IQIndexDat(i1,q1,i2,q2,i3,q3,i4,q4,i5,q5))
    {
    if(i2.type() != i1.type() 
    || i3.type() != i1.type()
    || i4.type() != i1.type()
    || i5.type() != i1.type())
        Error("Indices must have the same type");
    }

int
totalM(const IQIndexDat::StorageT& storage)
    {
    int tm = 0;
    Foreach(const IndexQN& iq, storage)
        {
        tm += iq.m();
#ifdef DEBUG
        if(iq.type() != storage.front().type())
            Error("Indices must have the same type");
#endif
        }
    return tm;
    }

IQIndex::
IQIndex(const string& name, 
        IQIndexDat::StorageT& ind_qn, 
        Arrow dir, int plev) 
    : 
    Index(name,totalM(ind_qn),ind_qn.front().type(),plev),
    _dir(dir), 
    pd(new IQIndexDat(ind_qn))
    { 
    }

IQIndex::
IQIndex(const IQIndex& other, 
        IQIndexDat::StorageT& ind_qn)
    : 
    Index(other.name(),totalM(ind_qn),other.type(),ind_qn.front().primeLevel()),
    _dir(other._dir), 
    pd(new IQIndexDat(ind_qn))
    { 
    }

IQIndex::
IQIndex(const Index& other, 
        const Index& i1, const QN& q1, 
        Arrow dir) 
    : 
    Index(other),
    _dir(dir), 
    pd(make_shared<IQIndexDat>(i1,q1))
    {
    Index::primeLevel(i1.primeLevel());
    }

IQIndex::
IQIndex(const Index& index, const IQIndexDatPtr& pdat)
    : 
    Index(index),
    _dir(In),
    pd(pdat)
    { }

void IQIndex::
write(ostream& s) const
    {
    IQINDEX_CHECK_NULL
    Index::write(s);
    s.write((char*)&_dir,sizeof(_dir));
    pd->write(s);
    }

void IQIndex::
read(istream& s)
    {
    Index::read(s);
    s.read((char*)&_dir,sizeof(_dir));
    pd = make_shared<IQIndexDat>();
    pd->read(s);
    }

string
showm(const IQIndex& I)
    {
#ifdef DEBUG
    if(I.isNull()) Error("Null IQIndex");
#endif
    string res = " ";
    ostringstream oh; 
    oh << I.m() << " | ";
    Foreach(const IndexQN& iq, I.indices())
        {
        oh << boost::format("[%d,%d,%s]:%d ") % iq.qn.sz() % iq.qn.Nf() % (iq.qn.sign()==1?"+":"-") % iq.m(); 
        }
    return oh.str();
    }

QN IQIndex::
qn(const Index& i) const
    { 
    IQINDEX_CHECK_NULL
    Foreach(const IndexQN& iq, *pd)
        { 
        if(iq == i) 
            return iq.qn; 
        }
    cerr << *this << "\n";
    cerr << "i = " << i << "\n";
    Error("IQIndex::qn(Index): IQIndex does not contain given index.");
    return QN();
    }

Arrow IQIndex::
dir() const { return _dir; }

void IQIndex::
conj() { _dir = -_dir; }

const Index& IQIndex::
findbyqn(QN q) const
    { 
    IQINDEX_CHECK_NULL
    Foreach(const IndexQN& iq, *pd)
        {
        if(iq.qn == q) 
            return iq;
        }
    Error("IQIndex::findbyqn: no Index had a matching QN.");
    return Index::Null();
    }

bool IQIndex::
hasindex(const Index& i) const
    { 
    IQINDEX_CHECK_NULL
    Foreach(const IndexQN& iq, *pd)
        {
        if(iq == i) 
            return true;
        }
    return false;
    }

bool IQIndex::
hasindex_noprime(const Index& i) const
    { 
    IQINDEX_CHECK_NULL
    Foreach(const IndexQN& iq, *pd)
        {
        if(iq.noprimeEquals(i)) 
            return true;
        }
    return false;
    }

void IQIndex::
primeLevel(int val)
    {
    solo();
    Index::primeLevel(val);
    Foreach(IndexQN& iq, *pd)
        iq.primeLevel(val);
    }

void IQIndex::
prime(int inc)
    {
    solo();
    Index::prime(inc);
    Foreach(IndexQN& iq, *pd)
        iq.prime(inc);
    }

void IQIndex::
prime(IndexType type, int inc)
    {
    solo();
    Index::prime(type,inc);
    Foreach(IndexQN& iq, *pd)
        iq.prime(type,inc);
    }

void IQIndex::
mapprime(int plevold, int plevnew, IndexType type)
    {
    solo();
    Index::mapprime(plevold,plevnew,type);
    Foreach(IndexQN& iq, *pd)
        iq.mapprime(plevold,plevnew,type);
    }

void IQIndex::
noprime(IndexType type)
    {
    solo();
    Index::noprime(type);
    Foreach(IndexQN& iq, *pd)
        iq.noprime(type);
    }


void IQIndex::
solo()
    {
    IQINDEX_CHECK_NULL
    if(!pd.unique())
        {
        const IQIndexDat& olddat = *pd;
        pd = make_shared<IQIndexDat>();
        pd->makeCopyOf(olddat);
        }
    }

const IQIndex& IQIndex::
Null()
    {
    static const IQIndex Null_(Index::Null(),IQIndexDat::Null());
    return Null_;
    }

const IQIndex& IQIndex::
IndReIm()
    {
    static const IQIndex IndReIm_(Index::IndReIm(),IQIndexDat::ReImDat());
    return IndReIm_;
    }

const IQIndex& IQIndex::
IndReImP()
    {
    static const IQIndex IndReImP_(Index::IndReImP(),IQIndexDat::ReImDatP());
    return IndReImP_;
    }

const IQIndex& IQIndex::
IndReImPP()
    {
    static const IQIndex IndReImPP_(Index::IndReImPP(),IQIndexDat::ReImDatPP());
    return IndReImPP_;
    }




IQIndexVal::
IQIndexVal()
    : iqind(IQIndex::Null()), i(1) 
    { }


IQIndexVal::
IQIndexVal(const IQIndex& iqindex, int i_) 
    : 
    iqind(iqindex),
    i(i_) 
    { 
    if(i > iqind.m() || i < 1) 
        {
        Print(iqindex);
        Print(i);
        Error("IQIndexVal: i out of range");
        }
    }


Index IQIndexVal::
index() const 
    { 
    int j,ii;
    calc_ind_ii(j,ii);
    return iqind.index(j);
    }


QN IQIndexVal::
qn() const 
    { 
    int j,ii;
    calc_ind_ii(j,ii);
    return iqind.qn(j);
    }

bool IQIndexVal::
operator==(const IQIndexVal& other) const
    {
    return (iqind == other.iqind && i == other.i);
    }

IQIndexVal::
operator IndexVal() const 
    { 
    return IndexVal(Index(iqind),i); 
    }


IndexVal IQIndexVal::
blockIndexVal() const 
    { 
    if(*this == IQIndexVal::Null())
        return IndexVal::Null();
    int j,ii;
    calc_ind_ii(j,ii);
    return IndexVal(iqind.index(j),ii); 
    }

/*

IQIndexVal::
operator ITensor() const 
    { 
    return ITensor(IndexVal(iqind,i)); 
    }
*/


void IQIndexVal::
calc_ind_ii(int& j, int& ii) const
    {
    j = 1;
    ii = i;
    while(ii > iqind.index(j).m())
        {
        ii -= iqind.index(j).m();
        ++j;
        }
    }

IQIndexVal IQIndex::
operator()(int n) const 
    { 
    return IQIndexVal(*this,n); 
    }

int
offset(const IQIndex& I, const Index& i)
    {
    int os = 0;
    Foreach(const IndexQN& iq, I.indices())
        {
        if(iq == i) return os;
        os += iq.m();
        }
    Print(I);
    Print(i);
    Error("Index not contained in IQIndex");
    return 0;
    }

ostream& 
operator<<(ostream &o, const IQIndex& I)
    {
    if(I.isNull()) 
        { 
        o << "IQIndex: (null)"; 
        return o;
        }
    o << "IQIndex: " << Index(I) << " <" << I.dir() << ">" << endl;
    for(int j = 1; j <= I.nindex(); ++j) 
        o << " " << I.index(j) SP I.qn(j) << "\n";
    return o;
    }

std::ostream& 
operator<<(std::ostream &s, const IndexQN& x)
    { 
    const Index& i = x;
    return s << "IndexQN: " << i
             << " (" << x.qn << ")\n";
    }

std::ostream& 
operator<<(std::ostream& s, const IQIndexVal& iv)
    { 
    return s << "IQIndexVal: i = " << iv.i 
             << ", iqind = " << iv.iqind << "\n"; 
    }
