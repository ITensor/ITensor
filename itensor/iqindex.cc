//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqindex.h"

using namespace std;
using boost::format;
using boost::array;

void 
intrusive_ptr_add_ref(IQIndexDat* pd) 
    { 
    ++(pd->numref); 
    }

void 
intrusive_ptr_release(IQIndexDat* pd) 
    { 
    if(!pd->is_static_ && --(pd->numref) == 0)
        { 
        delete pd; 
        } 
    }

IQIndexDat::
IQIndexDat() 
    : 
    numref(0), 
    is_static_(false) 
    { }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    iq_.push_back(inqn(i3,q3));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3,
           const Index& i4, const QN& q4)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    iq_.push_back(inqn(i3,q3));
    iq_.push_back(inqn(i4,q4));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3,
           const Index& i4, const QN& q4,
           const Index& i5, const QN& q5)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    iq_.push_back(inqn(i3,q3));
    iq_.push_back(inqn(i4,q4));
    iq_.push_back(inqn(i5,q5));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3,
           const Index& i4, const QN& q4,
           const Index& i5, const QN& q5,
           const Index& i6, const QN& q6)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    iq_.push_back(inqn(i3,q3));
    iq_.push_back(inqn(i4,q4));
    iq_.push_back(inqn(i5,q5));
    iq_.push_back(inqn(i6,q6));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3,
           const Index& i4, const QN& q4,
           const Index& i5, const QN& q5,
           const Index& i6, const QN& q6,
           const Index& i7, const QN& q7)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    iq_.push_back(inqn(i3,q3));
    iq_.push_back(inqn(i4,q4));
    iq_.push_back(inqn(i5,q5));
    iq_.push_back(inqn(i6,q6));
    iq_.push_back(inqn(i7,q7));
    }

IQIndexDat::
IQIndexDat(const Index& i1, const QN& q1,
           const Index& i2, const QN& q2,
           const Index& i3, const QN& q3,
           const Index& i4, const QN& q4,
           const Index& i5, const QN& q5,
           const Index& i6, const QN& q6,
           const Index& i7, const QN& q7,
           const Index& i8, const QN& q8)
    : 
    numref(0), 
    is_static_(false)
    {
    iq_.push_back(inqn(i1,q1));
    iq_.push_back(inqn(i2,q2));
    iq_.push_back(inqn(i3,q3));
    iq_.push_back(inqn(i4,q4));
    iq_.push_back(inqn(i5,q5));
    iq_.push_back(inqn(i6,q6));
    iq_.push_back(inqn(i7,q7));
    iq_.push_back(inqn(i8,q8));
    }

IQIndexDat::
IQIndexDat(std::vector<inqn>& ind_qn)
    : 
    numref(0), 
    is_static_(false)
    { 
    iq_.swap(ind_qn); 
    }

IQIndexDat::
IQIndexDat(const IQIndexDat& other) 
    : 
    iq_(other.iq_),
    numref(0), 
    is_static_(other.is_static_)
    { 
    }

IQIndexDat::
IQIndexDat(std::istream& s) 
    : 
    numref(0), 
    is_static_(false) 
    { read(s); }

IQIndexDat::
IQIndexDat(Index::Imaker im)
    : 
    numref(0), 
    is_static_(true)
    { 
    if(im == Index::makeNull)
        iq_.push_back(inqn(Index::Null(),QN())); 
    else if(im == Index::makeReIm)
        iq_.push_back(inqn(Index::IndReIm(),QN())); 
    else if(im == Index::makeReImP)
        iq_.push_back(inqn(Index::IndReImP(),QN())); 
    else if(im == Index::makeReImPP)
        iq_.push_back(inqn(Index::IndReImPP(),QN())); 
    }

void IQIndexDat::
write(std::ostream& s) const
    {
    size_t size = iq_.size();
    s.write((char*)&size,sizeof(size));
    for(const_iq_it x = iq_.begin(); x != iq_.end(); ++x)
        { x->write(s); }
    }

void IQIndexDat::
read(std::istream& s)
    {
    size_t size; s.read((char*)&size,sizeof(size));
    iq_.resize(size);
    for(iq_it x = iq_.begin(); x != iq_.end(); ++x)
        { x->read(s); }
    numref = 0;
    }

IQIndexDat* IQIndexDat::
Null()
    {
    static IQIndexDat Null_(Index::makeNull);
    return &Null_;
    }

IQIndexDat* IQIndexDat::
ReImDat()
    {
    static IQIndexDat ReImDat_(Index::makeReIm);
    return &ReImDat_;
    }

IQIndexDat* IQIndexDat::
ReImDatP()
    {
    static IQIndexDat ReImDatP_(Index::makeReImP);
    return &ReImDatP_;
    }

IQIndexDat* IQIndexDat::
ReImDatPP()
    {
    static IQIndexDat ReImDatPP_(Index::makeReImPP);
    return &ReImDatPP_;
    }

//
// IQIndex Methods
//

#ifdef DEBUG
#define IQINDEX_CHECK_NULL if(pd == 0) Error("IQIndex is null");
#else
#define IQINDEX_CHECK_NULL
#endif

const std::vector<inqn>& IQIndex::
iq() const 
    { 
    IQINDEX_CHECK_NULL
    return pd->iq_;
    }

int IQIndex::
nindex() const 
    { 
    IQINDEX_CHECK_NULL
    return (int) pd->iq_.size(); 
    }

const Index& IQIndex::
index(int i) const 
    {
    IQINDEX_CHECK_NULL
    return GET(pd->iq_,i-1).index; 
    }

const QN& IQIndex::
qn(int i) const 
    {
    IQINDEX_CHECK_NULL
    return GET(pd->iq_,i-1).qn; 
    }

IQIndex::
IQIndex() 
    : 
    _dir(Out), 
    pd(0) 
    { }

IQIndex::
IQIndex(const Index& other, Arrow dir)
    : 
    index_(other), 
    _dir(dir), 
    pd(0)
    { }

IQIndex::
IQIndex(const std::string& name,
                 IndexType it, 
                 Arrow dir, 
                 int plev) 
    : 
    index_(name,1,it,plev), 
    _dir(dir), 
    pd(0) 
    { }

IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        Arrow dir) 
    : 
    index_(name,i1.m(),i1.type(),i1.primeLevel()), 
    _dir(dir), 
    pd(new IQIndexDat(i1,q1))
    {
    }

IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        Arrow dir) 
    : 
    index_(name,i1.m()+i2.m(),i1.type(),i1.primeLevel()), 
    _dir(dir), 
    pd(new IQIndexDat(i1,q1,i2,q2))
    {
    if(i2.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        const Index& i3, const QN& q3,
        Arrow dir) 
    : 
    index_(name,i1.m()+i2.m()+i3.m(),i1.type(),i1.primeLevel()), 
    _dir(dir),
    pd(new IQIndexDat(i1,q1,i2,q2,i3,q3))
    {
    if(i2.type() != i1.type() 
    || i3.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        const Index& i3, const QN& q3,
        const Index& i4, const QN& q4,
        Arrow dir) 
    : 
    index_(name,i1.m()+i2.m()+i3.m()+i4.m(),i1.type(),i1.primeLevel()), 
    _dir(dir),
    pd(new IQIndexDat(i1,q1,i2,q2,i3,q3,i4,q4))
    {
    if(i2.type() != i1.type() 
    || i3.type() != i1.type()
    || i4.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const std::string& name, 
        const Index& i1, const QN& q1, 
        const Index& i2, const QN& q2,
        const Index& i3, const QN& q3,
        const Index& i4, const QN& q4,
        const Index& i5, const QN& q5,
        Arrow dir) 
    : 
    index_(name,i1.m()+i2.m()+i3.m()+i4.m()+i5.m(),i1.type(),i1.primeLevel()), 
    _dir(dir),
    pd(new IQIndexDat(i1,q1,i2,q2,i3,q3,i4,q4,i5,q5))
    {
    if(i2.type() != i1.type() 
    || i3.type() != i1.type()
    || i4.type() != i1.type()
    || i5.type() != i1.type())
        Error("Indices must have the same type");
    }

IQIndex::
IQIndex(const std::string& name, 
        std::vector<inqn>& ind_qn, 
        Arrow dir, int plev) 
    : 
    _dir(dir), 
    pd(new IQIndexDat(ind_qn))
    { 
    int mm = 0;
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        {
        mm += x->index.m();
        if(x->index.type() != index(1).type())
            Error("Indices must have the same type");
        }
    index_ = Index(name,mm,index(1).type(),plev);
    }

IQIndex::
IQIndex(const IQIndex& other, 
        std::vector<inqn>& ind_qn)
    : 
    _dir(other._dir), 
    pd(new IQIndexDat(ind_qn))
    { 
    int mm = 0;
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        {
        mm += x->index.m();
        if(x->index.type() != index(1).type())
            Error("Indices must have the same type");
        }
    index_ = Index(other.name(),mm,other.type(),
                   pd->iq_.back().index.primeLevel()); 
    }

IQIndex::
IQIndex(const Index& other, 
        const Index& i1, const QN& q1, 
        Arrow dir) 
    : 
    index_(other),
    _dir(dir), 
    pd(new IQIndexDat(i1,q1))
    {
    index_.primeLevel(i1.primeLevel());
    }

IQIndex::
IQIndex(IndexType type, const IQIndex& other, int inc) 
    : 
    index_(other), 
    _dir(other._dir), 
    pd(other.pd)  
    { doprime(type,inc); }

IQIndex::
IQIndex(std::istream& s) 
    { read(s); }

IQIndex::
IQIndex(Index::Imaker im)
    : 
    index_(im), 
    _dir(In)
    {
    if(im == Index::makeNull)
        { pd = IQIndexDat::Null(); }
    else if(im == Index::makeReIm)
        { pd = IQIndexDat::ReImDat(); }
    else if(im == Index::makeReImP)
        { pd = IQIndexDat::ReImDatP(); }
    else if(im == Index::makeReImPP)
        { pd = IQIndexDat::ReImDatPP(); }
    else Error("IQIndex: Unrecognized Imaker type.");
    }

void IQIndex::
write(std::ostream& s) const
    {
    IQINDEX_CHECK_NULL
    index_.write(s);
    s.write((char*)&_dir,sizeof(_dir));
    pd->write(s);
    }

void IQIndex::
read(std::istream& s)
    {
    index_.read(s);
    s.read((char*)&_dir,sizeof(_dir));
    pd = new IQIndexDat(s);
    }

int IQIndex::
biggestm() const
    {
    IQINDEX_CHECK_NULL
    int mm = 0;
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        { mm = max(mm,x->index.m()); }
    return mm;
    }

std::string IQIndex::
showm() const
    {
    IQINDEX_CHECK_NULL
    std::string res = " ";
    std::ostringstream oh; 
    oh << this->m() << " | ";
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        {
        QN q = x->qn;
        oh << boost::format("[%d,%d,%s]:%d ") % q.sz() % q.Nf() % (q.sign()==1?"+":"-") % x->index.m(); 
        }
    return oh.str();
    }

void IQIndex::
negate()
    {
    IQINDEX_CHECK_NULL
    for(iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        { x->qn *= -1; }
    }

IQIndex 
negate(IQIndex I) // Quantum numbers negated
    { 
    I.negate();
    return I;
    }

QN IQIndex::
qn(const Index& i) const
    { 
    IQINDEX_CHECK_NULL
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        { if(x->index == i) return x->qn; }
    std::cerr << *this << "\n";
    std::cerr << "i = " << i << "\n";
    Error("IQIndex::qn(Index): IQIndex does not contain given index.");
    return QN();
    }

Arrow IQIndex::
dir() const { return _dir; }

void IQIndex::
conj() { _dir = _dir*Switch; }

const Index& IQIndex::
findbyqn(QN q) const
    { 
    IQINDEX_CHECK_NULL
    for(size_t i = 0; i < pd->iq_.size(); ++i)
        if(pd->iq_[i].qn == q) return pd->iq_[i].index;
    Error("IQIndex::findbyqn: no Index had a matching QN.");
    return pd->iq_[0].index;
    }

bool IQIndex::
hasindex(const Index& i) const
    { 
    IQINDEX_CHECK_NULL
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        if(x->index == i) return true;
    return false;
    }

bool IQIndex::
hasindex_noprime(const Index& i) const
    { 
    IQINDEX_CHECK_NULL
    for(const_iq_it x = pd->iq_.begin(); x != pd->iq_.end(); ++x)
        if(x->index.noprime_equals(i)) return true;
    return false;
    }

int IQIndex::
offset(const Index& I) const
    {
    int os = 0;
    for(size_t j = 0; j < pd->iq_.size(); ++j)
        {
        const Index& J = pd->iq_[j].index;
        if(J == I) return os;
        os += J.m();
        }
    Print(*this);
    Print(I);
    Error("Index not contained in IQIndex");
    return 0;
    }

void IQIndex::
doprime(IndexType type, int inc)
    {
    solo();
    index_.doprime(type,inc);
    DoPrimer dp(type,inc);
    for_each(pd->iq_.begin(),pd->iq_.end(),dp);
    }

void IQIndex::
mapprime(int plevold, int plevnew, IndexType type)
    {
    solo();
    index_.mapprime(plevold,plevnew,type);
    for_each(pd->iq_.begin(),pd->iq_.end(),MapPrimer(plevold,plevnew,type));
    }

void IQIndex::
noprime(IndexType type)
    {
    solo();
    index_.noprime(type);
    for(size_t j = 0; j < pd->iq_.size(); ++j)
    { pd->iq_[j].index.noprime(type); }
    }

IQIndex IQIndex::
primed(int inc) const
    {
    return IQIndex(All,*this,inc);
    }

void IQIndex::
print(std::string name) const
    { std::cerr << "\n" << name << " =\n" << *this << "\n"; }

void IQIndex::
solo()
    {
    IQINDEX_CHECK_NULL
    if(pd->count() != 1)
        {
        pd = new IQIndexDat(*pd);
        }
    }

const IQIndex& IQIndex::
Null()
    {
    static const IQIndex Null_(Index::makeNull);
    return Null_;
    }

const IQIndex& IQIndex::
IndReIm()
    {
    static const IQIndex IndReIm_(Index::makeReIm);
    return IndReIm_;
    }

const IQIndex& IQIndex::
IndReImP()
    {
    static const IQIndex IndReImP_(Index::makeReImP);
    return IndReImP_;
    }

const IQIndex& IQIndex::
IndReImPP()
    {
    static const IQIndex IndReImPP_(Index::makeReImPP);
    return IndReImPP_;
    }


std::ostream& 
operator<<(std::ostream &o, const IQIndex& I)
    {
    if(I.isNull()) 
        { 
        o << "IQIndex: (null)"; 
        return o;
        }
    o << "IQIndex: " << I.index_ << " <" << I.dir() << ">" << std::endl;
    for(int j = 1; j <= I.nindex(); ++j) 
        o << " " << I.index(j) SP I.qn(j) << "\n";
    return o;
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
#ifdef DEBUG
    if(iqindex == IQIndex::Null())
        {
        Error("IQIndexVal index set to null");
        }
#endif
    if(i > iqind.m() || i < 1) 
        {
        Print(iqindex);
        Print(i);
        Error("IQIndexVal: i out of range");
        }
    }

IQIndexVal::
IQIndexVal(Index::Imaker im)
    {
    if(im == Index::makeNull)
        {
        iqind = IQIndex::Null();
        i = 1;
        }
    else
        {
        Error("Only makeNull supported in Imaker constructor");
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
    { return IQIndexVal(*this,n); }
