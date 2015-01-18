//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqindex.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ostringstream;

//
// IQIndexDat
//


void IQIndexDat::
write(ostream& s) const
    {
    size_t size = iq_.size();
    s.write((char*)&size,sizeof(size));
    for(const IndexQN& x : iq_)
        { 
        x.write(s); 
        }
    }

void IQIndexDat::
read(istream& s)
    {
    size_t size; s.read((char*)&size,sizeof(size));
    iq_.resize(size);
    for(IndexQN& x : iq_)
        { 
        x.read(s); 
        }
    }

const IQIndexDatPtr& IQIndexDat::
Null()
    {
    static IQIndexDatPtr Null_ = make_shared<IQIndexDat>(Index(),QN());
    return Null_;
    }

//
// IQIndex Methods
//

#ifdef DEBUG
#define IQINDEX_CHECK_NULL if(pd == 0) throw ITError("IQIndex is null");
#else
#define IQINDEX_CHECK_NULL
#endif

const IQIndexDat::storage& IQIndex::
indices() const 
    { 
    IQINDEX_CHECK_NULL
    return pd->indices();
    }

long IQIndex::
nindex() const 
    { 
    IQINDEX_CHECK_NULL
    return (long) pd->size(); 
    }

const Index& IQIndex::
index(long i) const 
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

const Index& IQIndex::
operator[](long i) const 
    {
    IQINDEX_CHECK_NULL
#ifdef DEBUG
    if(i > nindex()-1)
        {
        Print(nindex());
        Print(i);
        Error("IQIndex::index arg out of range");
        }
#endif
    return pd->operator[](i);
    }

const QN& IQIndex::
qn(long i) const 
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

long
totalM(const IQIndexDat::storage& storage)
    {
    long tm = 0;
    for(const IndexQN& iq : storage)
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
IQIndex(const Index& index, const IQIndexDatPtr& pdat)
    : 
    Index(index),
    pd(pdat),
    dir_(In)
    { }

IQIndex& IQIndex::
dag() 
    { 
    dir_ = -dir_; 
    return *this;
    }

void IQIndex::
write(ostream& s) const
    {
    IQINDEX_CHECK_NULL
    Index::write(s);
    s.write((char*)&dir_,sizeof(dir_));
    pd->write(s);
    }

IQIndex& IQIndex::
read(istream& s)
    {
    Index::read(s);
    s.read((char*)&dir_,sizeof(dir_));
    pd = make_shared<IQIndexDat>();
    pd->read(s);
    return *this;
    }

string
showm(const IQIndex& I)
    {
#ifdef DEBUG
    if(!I) Error("Calling showm on null IQIndex");
#endif
    string res = " ";
    ostringstream oh; 
    oh << I.m() << " | ";
    for(const IndexQN& iq : I.indices())
        {
        oh << iq.qn << ":" << iq.m() << " ";
        }
    return oh.str();
    }


IQIndex& IQIndex::
primeLevel(int val)
    {
    solo();
    Index::primeLevel(val);
    for(IndexQN& iq : *pd)
        {
        iq.primeLevel(val);
        }
    return *this;
    }

IQIndex& IQIndex::
prime(int inc)
    {
    solo();
    Index::prime(inc);
    for(IndexQN& iq : *pd)
        iq.prime(inc);
    return *this;
    }

IQIndex& IQIndex::
prime(IndexType type, int inc)
    {
    solo();
    Index::prime(type,inc);
    for(IndexQN& iq : *pd)
        iq.prime(type,inc);
    return *this;
    }

IQIndex& IQIndex::
mapprime(int plevold, int plevnew, IndexType type)
    {
    solo();
    Index::mapprime(plevold,plevnew,type);
    for(IndexQN& iq : *pd)
        iq.mapprime(plevold,plevnew,type);
    return *this;
    }

IQIndex& IQIndex::
noprime(IndexType type)
    {
    solo();
    Index::noprime(type);
    for(IndexQN& iq : *pd)
        iq.noprime(type);
    return *this;
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

//const IQIndex& IQIndex::
//Null()
//    {
//    static const IQIndex Null_(Index::Null(),IQIndexDat::Null());
//    return Null_;
//    }

void
calc_ind_ii(const IQIndexVal& iv, long& j, long& ii)
    {
    j = 1;
    ii = iv.i;
    while(ii > iv.index.index(j).m())
        {
        ii -= iv.index.index(j).m();
        ++j;
        }
    }


IQIndexVal::
IQIndexVal()
    : 
    i(0) 
    { }


IQIndexVal::
IQIndexVal(const IQIndex& iqindex, long i_) 
    : 
    index(iqindex),
    i(i_) 
    { 
#ifdef DEBUG
    if(i > m() || i < 1) 
        {
        Print(index);
        Print(i);
        Error("IQIndexVal: i out of range");
        }
#endif
    }


IndexQN IQIndexVal::
indexqn() const 
    { 
    long j,ii;
    calc_ind_ii(*this,j,ii);
    return IndexQN(index.index(j),index.qn(j));
    }


const QN& IQIndexVal::
qn() const 
    { 
    long j,ii;
    calc_ind_ii(*this,j,ii);
    return index.qn(j);
    }

bool IQIndexVal::
operator==(const IQIndexVal& other) const
    {
    return (index == other.index && i == other.i);
    }

IQIndexVal::
operator IndexVal() const 
    { 
    return IndexVal(Index(index),i); 
    }


IndexVal IQIndexVal::
blockIndexVal() const 
    { 
    //if(*this == IQIndexVal::Null())
    //    return IndexVal::Null();
    long j,ii;
    calc_ind_ii(*this,j,ii);
    return IndexVal(index.index(j),ii); 
    }

IQIndexVal&  IQIndexVal::
prime(int inc)
    {
    index.prime(inc);
    return *this;
    }

IQIndexVal&  IQIndexVal::
prime(IndexType type, int inc)
    {
    index.prime(type,inc);
    return *this;
    }

IQIndexVal&  IQIndexVal::
noprime(IndexType type)
    {
    index.noprime(type);
    return *this;
    }

IQIndexVal&  IQIndexVal::
mapprime(int plevold, int plevnew, IndexType type)
    {
    index.mapprime(plevold,plevnew,type);
    return *this;
    }

IQIndexVal&  IQIndexVal::
dag() { index.dag(); return *this; }

/*

IQIndexVal::
operator ITensor() const 
    { 
    return ITensor(IndexVal(iqind,i)); 
    }
*/



IQIndexVal IQIndex::
operator()(long n) const 
    { 
    return IQIndexVal(*this,n); 
    }

bool
hasindex(const IQIndex& J, const Index& i)
    { 
    for(const Index& j : J.indices())
        {
        if(j == i) return true;
        }
    return false;
    }

long
findindex(const IQIndex& J, const Index& i)
    { 
    for(long j = 1; j <= J.nindex(); ++j)
        {
        if(J.index(j) == i) return j;
        }
    return 0;
    }

long
offset(const IQIndex& I, const Index& i)
    {
    long os = 0;
    for(const IndexQN& iq : I.indices())
        {
        if(iq == i) return os;
        os += iq.m();
        }
    Print(I);
    Print(i);
    Error("Index not contained in IQIndex");
    return 0;
    }

QN
qn(const IQIndex& I, const Index& i)
    { 
    for(const IndexQN& jq : I.indices())
        { 
        if(jq == i) 
            return jq.qn; 
        }
    cout << I << "\n";
    cout << "i = " << i << endl;
    Error("IQIndex does not contain given index.");
    return QN();
    }

Index
findByQN(const IQIndex& I, const QN& qn)
    { 
    for(const IndexQN& jq : I.indices())
        { 
        if(jq.qn == qn) 
            return jq;
        }
    cout << I << "\n";
    cout << "qn = " << qn << endl;
    Error("IQIndex does not contain given QN block.");
    return Index();
    }

ostream& 
operator<<(ostream &o, const IQIndex& I)
    {
    if(!I) 
        { 
        o << "IQIndex: (null)"; 
        return o;
        }
    o << "IQIndex: " << Index(I) << " <" << I.dir() << ">" << endl;
    for(long j = 1; j <= I.nindex(); ++j) 
        o << "  " << I.index(j) SP I.qn(j) << "\n";
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
    const IQIndex& I = iv.index;
    return s << "IQIndexVal: i = " << iv.i << " for IQIndex:\n  " << I << "\n"; 
    }

}; //namespace itensor
