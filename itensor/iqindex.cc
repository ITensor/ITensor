//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itensor.h"
#include "itensor/iqindex.h"
#include "itensor/util/print_macro.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::vector;
using std::string;
using std::ostringstream;
using std::make_shared;

//
// IQIndex Methods
//

#ifdef DEBUG
#define IQINDEX_CHECK_NULL if(pd == 0) throw ITError("IQIndex storage unallocated");
#else
#define IQINDEX_CHECK_NULL
#endif

//const IQIndexDat::storage& IQIndex::
//inds() const 
//    { 
//    IQINDEX_CHECK_NULL
//    return pd->inds();
//    }

IQIndex::const_iterator IQIndex::
begin() const 
    { 
    IQINDEX_CHECK_NULL
    return pd->begin();
    }

IQIndex::const_iterator IQIndex::
end() const 
    { 
    IQINDEX_CHECK_NULL
    return pd->end();
    }

long IQIndex::
nindex() const 
    { 
    IQINDEX_CHECK_NULL
    return (long) pd->size(); 
    }

Index IQIndex::
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
    return itensor::prime(pd->index(i),Index::primeLevel());
    }

Index IQIndex::
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
    return itensor::prime(pd->operator[](i),Index::primeLevel());
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
        tm += iq.index.m();
#ifdef DEBUG
        if(iq.index.type() != storage.front().type())
            Error("Indices must have the same type");
#endif
        }
    return tm;
    }

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
    itensor::write(s,dir_);
    itensor::write(s,*pd);
    }

IQIndex& IQIndex::
read(istream& s)
    {
    Index::read(s);
    itensor::read(s,dir_);
    pd = make_shared<IQIndexDat>();
    itensor::read(s,*pd);
    return *this;
    }

string
showm(IQIndex const& I)
    {
#ifdef DEBUG
    if(!I) Error("Calling showm on null IQIndex");
#endif
    ostringstream oh; 
    oh << "m=" << I.m() << " | ";
    for(const IndexQN& iq : I)
        {
        oh << iq.qn << ":" << iq.m() << " ";
        }
    return oh.str();
    }


//IQIndex& IQIndex::
//primeLevel(int val)
//    {
//    //solo();
//    Index::primeLevel(val);
//    //for(IndexQN& iq : *pd)
//    //    {
//    //    iq.primeLevel(val);
//    //    }
//    return *this;
//    }

//IQIndex& IQIndex::
//prime(int inc)
//    {
//    //solo();
//    Index::prime(inc);
//    //for(IndexQN& iq : *pd)
//    //    iq.prime(inc);
//    return *this;
//    }
//
//IQIndex& IQIndex::
//prime(IndexType type, int inc)
//    {
//    //solo();
//    Index::prime(type,inc);
//    //for(IndexQN& iq : *pd)
//    //    iq.prime(type,inc);
//    return *this;
//    }
//
//IQIndex& IQIndex::
//mapprime(int plevold, int plevnew, IndexType type)
//    {
//    //solo();
//    Index::mapprime(plevold,plevnew,type);
//    //for(IndexQN& iq : *pd)
//    //    iq.mapprime(plevold,plevnew,type);
//    return *this;
//    }
//
//IQIndex& IQIndex::
//noprime(IndexType type)
//    {
//    //solo();
//    Index::noprime(type);
//    //for(IndexQN& iq : *pd)
//    //    iq.noprime(type);
//    return *this;
//    }


//void IQIndex::
//solo()
//    {
//    IQINDEX_CHECK_NULL
//    if(!pd.unique()) pd = pd->clone();
//    }

void
calc_ind_ii(const IQIndexVal& iv, long& j, long& ii)
    {
    j = 1;
    ii = iv.val;
    while(ii > iv.index.index(j).m())
        {
        ii -= iv.index.index(j).m();
        ++j;
        }
    }


IQIndexVal::
IQIndexVal()
    : 
    val(0) 
    { }


IQIndexVal::
IQIndexVal(const IQIndex& iqindex, long val_) 
    : 
    index(iqindex),
    val(val_) 
    { 
#ifdef DEBUG
    //if(val > m() || val < 1) 
    //    {
    //    Print(index);
    //    Print(val);
    //    Error("IQIndexVal: val out of range");
    //    }
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

bool
operator==(IQIndexVal const& iv1, IQIndexVal const& iv2)
    {
    return (iv1.index == iv2.index && iv1.val == iv2.val);
    }

IQIndexVal::
operator IndexVal() const 
    { 
    return IndexVal(Index(index),val); 
    }


IndexVal IQIndexVal::
blockIndexVal() const 
    { 
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

ITensor
operator*(IQIndexVal const& iqiv, IndexVal const& iv)
    { 
    return IndexVal(Index(iqiv.index),iqiv.val) * iv; 
    }

/*

IQIndexVal::
operator ITensor() const 
    { 
    return ITensor(IndexVal(iqind,i)); 
    }
*/

IQIndexVal IQIndex::
operator()(long val) const 
    { 
    return IQIndexVal(*this,val); 
    }

bool
hasindex(IQIndex const& J, Index const& i)
    { 
    for(auto& j : J)
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
    for(const IndexQN& iq : I)
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
    for(const IndexQN& jq : I)
        { 
        if(jq == i) 
            return jq.qn; 
        }
    println("I = ",I);
    println("i = ",i);
    Error("IQIndex does not contain given index.");
    return QN();
    }

Index
findByQN(const IQIndex& I, const QN& qn)
    { 
    for(auto& iq : I)
        { 
        if(iq.qn == qn) return Index(iq);
        }
    println("I = ",I);
    println("qn = ",qn);
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
    o << "IQIndex" << Index(I) << " <" << I.dir() << ">" << "\n";
    for(long j = 1; j <= I.nindex(); ++j) 
        o << "  " << I.index(j) << " " <<  I.qn(j) << "\n";
    return o;
    }

std::ostream& 
operator<<(std::ostream &s, IndexQN const& x)
    { 
    return s << "IndexQN: " << x.index
             << " (" << x.qn << ")\n";
    }

std::ostream& 
operator<<(std::ostream& s, IQIndexVal const& iv)
    { 
    return s << "IQIndexVal: val = " << iv.val << " for IQIndex:\n  " << iv.index << "\n"; 
    }

} //namespace itensor
