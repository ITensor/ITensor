//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/itensor.h"
#include "itensor/iqindex.h"
#include "itensor/util/print_macro.h"
#include "itensor/util/readwrite.h"

namespace itensor {

using std::istream;
using std::ostream;
using std::vector;
using std::string;
using std::ostringstream;
using std::make_shared;

long
totalM(IQIndex::storage const& storage)
    {
    long tm = 0;
    for(auto& iq : storage)
        {
        tm += iq.index.m();
#ifdef DEBUG
        if(iq.index.type() != storage.front().type())
            {
            Print(iq.index.type());
            Print(storage.front().type());
            Error("Indices must have the same type");
            }
#endif
        }
    return tm;
    }

class IQIndexDat
    {
    public:
    using storage = std::vector<IndexQN>;
    using iterator = storage::iterator;
    using const_iterator = storage::const_iterator;
    private:
    storage iq_;
    public:

    IQIndexDat() { }

    //template<typename... Rest>
    //IQIndexDat(Index const& i1, 
    //           QN const& q1,
    //           Rest const&... args) 
    //    { 
    //    constexpr auto size = sizeof...(args)/2+1;
    //    iq_ = stdx::reserve_vector<IndexQN>(size);
    //    detail::fill(iq_,i1,q1,args...);
    //    }

    explicit
    IQIndexDat(storage const& ind_qn) 
      : iq_(ind_qn)
        { }

    explicit
    IQIndexDat(storage&& ind_qn) 
      : iq_(std::move(ind_qn))
        { }

    //Disallow copying
    IQIndexDat(IQIndexDat const&) = delete;

    void 
    operator=(IQIndexDat const&) = delete;

    void
    setStore(storage && iq) { iq_ = std::move(iq); }

    storage const&
    inds() const { return iq_; }

    long
    size() const { return iq_.size(); }

    Index const&
    index(long i) { return iq_[i-1].index; }

    Index const&
    operator[](long i) { return iq_[i].index; }

    QN const&
    qn(long i) { return iq_[i-1].qn; }

    iterator
    begin() { return iq_.begin(); }

    iterator
    end() { return iq_.end(); }

    const_iterator
    begin() const { return iq_.begin(); }

    const_iterator
    end()   const { return iq_.end(); }

    storage const&
    store() const { return iq_; }
    };

void 
write(std::ostream & s, IQIndexDat const& d) 
    { 
    write(s,d.store()); 
    }

void 
read(std::istream & s, IQIndexDat & d) 
    { 
    IQIndexDat::storage store;
    read(s,store); 
    d.setStore(std::move(store));
    }

//
// IQIndex Methods
//

#ifdef DEBUG
#define IQINDEX_CHECK_NULL if(pd == 0) throw ITError("IQIndex storage unallocated");
#else
#define IQINDEX_CHECK_NULL
#endif

IQIndex::
IQIndex(std::string const& name, 
        storage && ind_qn, 
        Arrow dir, 
        int plev) 
  : Index(name,totalM(ind_qn),ind_qn.front().index.type(),plev),
    dir_(dir)
    { 
    makeStorage(std::move(ind_qn));
    }

// Constructor taking a storage pointer
IQIndex::
IQIndex(storage_ptr const& p,
        std::string const& name, 
        Arrow dir,
        int plev)
  : Index(name,totalM(p->store()),p->store().front().index.type(),plev),
    pd(p),
    dir_(dir)
    {
    }

//const IQIndexDat::storage& IQIndex::
//inds() const 
//    { 
//    IQINDEX_CHECK_NULL
//    return pd->inds();
//    }

//IQIndex::const_iterator IQIndex::
//begin() const 
//    { 
//    IQINDEX_CHECK_NULL
//    return pd->begin();
//    }
//
//IQIndex::const_iterator IQIndex::
//end() const 
//    { 
//    IQINDEX_CHECK_NULL
//    return pd->end();
//    }

long IQIndex::
nblock() const 
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

IQIndex
sim(IQIndex const& I, int plev)
    {
    return IQIndex(I.store(),"~"+I.rawname(),I.dir(),plev);
    }

string
showm(IQIndex const& I)
    {
#ifdef DEBUG
    if(!I) Error("Calling showm on null IQIndex");
#endif
    ostringstream oh; 
    oh << "m=" << I.m() << " | ";
    for(auto iq : I)
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

struct IndSector
    {
    long sector = 0l;
    long sind   = 0l;

    IndSector(long sec, long si) : sector(sec), sind(si) { }
    };

IndSector
sectorInfo(IQIndexVal const& iv)
    {
    auto is = IndSector(1,iv.val);
    while(is.sind > iv.index.index(is.sector).m())
        {
        is.sind -= iv.index.index(is.sector).m();
        is.sector += 1;
        }
    return is;
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
    auto is = sectorInfo(*this);
    return IndexQN(index.index(is.sector),index.qn(is.sector));
    }


const QN& IQIndexVal::
qn() const 
    { 
    auto is = sectorInfo(*this);
    return index.qn(is.sector);
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
    auto is = sectorInfo(*this);
    return IndexVal(index.index(is.sector),is.sind); 
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

void IQIndex::
makeStorage(storage && iq)
    {
    pd = std::make_shared<IQIndexDat>(std::move(iq));
    }

bool
hasindex(IQIndex const& J, Index const& i)
    { 
    for(auto n : range(J.nindex()))
        {
        if(J[n] == i) return true;
        }
    return false;
    }

long
findindex(IQIndex const& J, Index const& i)
    { 
    for(auto j : range1(J.nindex()))
        {
        if(J.index(j) == i) return j;
        }
    return 0;
    }

long
offset(IQIndex const& I, Index const& i)
    {
    long os = 0;
    for(auto n : range(I.nindex()))
        {
        if(I[n] == i) return os;
        os += I[n].m();
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
findByQN(IQIndex const& I, QN const& qn)
    { 
    for(auto iq : I)
        { 
        if(iq.qn == qn) return Index(iq);
        }
    println("I = ",I);
    println("qn = ",qn);
    Error("IQIndex does not contain given QN block.");
    return Index();
    }

ostream& 
operator<<(ostream & o, IQIndex const& I)
    {
    if(!I) 
        { 
        o << "IQIndex: (null)"; 
        return o;
        }
    o << "IQIndex" << Index(I) << " <" << I.dir() << ">" << "\n";
    for(auto j : range1(I.nindex()))
        {
        o << "  " << I.index(j) << " " <<  I.qn(j) << "\n";
        }
    return o;
    }

void IndexQN::
write(std::ostream & s) const
    { 
    index.write(s); 
    itensor::write(s,qn); 
    }

void IndexQN::
read(std::istream& s)
    { 
    index.read(s); 
    itensor::read(s,qn); 
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
