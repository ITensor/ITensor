//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/index.h"
#include "itensor/util/readwrite.h"
#include "itensor/util/print_macro.h"

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

// nameint() is a weird name without Index names, depecrate
//string 
//nameint(string const& f, int n)
//    { 
//    return format("%s%d",f,n);
//    }

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
Index(long m, 
      TagSet const& t)
  : id_(generateID()),
    m_(m),
    primelevel_(0),
    tags_(t)
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


IndexVal Index::
operator()(long val) const
    {
    return IndexVal(*this,val);
    }
IndexVal Index::
operator=(long val) const 
    { 
    return operator()(val); 
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
equalsIgnorePrime(Index const& i1, Index const& i2)
    { 
    return (i1.id() == i2.id()) && (tags(i1) == tags(i2));
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
operator<<(std::ostream & s, Index const& I)
    {
    s << "(";
    if(size(tags(I)) > 0) s << tags(I) << ",";
    s << I.m();
    if(Global::showIDs()) 
        {
        s << "|id=" << (I.id() % 1000);
        //s << "," << I.id();
        }
    s << ")"; 
    if(I.primeLevel() > 0) 
        {
        if(I.primeLevel() > 3)
            {
            s << "'" << I.primeLevel();
            }
        else
            {
            for(int n = 1; n <= I.primeLevel(); ++n)
                s << "'";
            }
        }
    if(hasQNs(I))
        {
        s << " <" << I.dir() << ">\n";
        for(auto j : range1(I.nblock()))
            {
            s << "  " << j << ": " << I.blocksize(j) << " " <<  I.qn(j);
            if(j != I.nblock()) s << "\n";
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
IndexVal(Index const& index_, 
         long val_) 
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
    return Index(I.m(),I.tags()).prime(plev);
    }

string
showm(Index const& I) { return format("m=%d",I.m()); }

std::ostream& 
operator<<(std::ostream& s, IndexVal const& iv)
    { 
    return s << "IndexVal: val = " << iv.val 
             << ", ind = " << iv.index;
    }

void
add(Args            & args, 
    Args::Name const& name, 
    TagSet     const& ts) 
    { 
    args.add(name,toString(ts)); 
    }

TagSet
getTagSet(Args       const& args, 
          Args::Name const& name)
    {
    if(!args.defined(name)) Error(format("Name %s not found in Args",name));
    return TagSet(args.getString(name));
    }

TagSet
getTagSet(const Args& args, 
          const Args::Name& name, 
          TagSet const& default_val)
    {
    if(!args.defined(name)) return default_val; 
    return TagSet(args.getString(name));
    }

struct IndSector
    {
    long sector = 0l;
    long sind   = 0l;

    IndSector(long sec, long si) : sector(sec), sind(si) { }
    };

IndSector
sectorInfo(IndexVal const& iv)
    {
    auto is = IndSector(1,iv.val);
    while(is.sind > iv.index.blocksize(is.sector))
        {
        is.sind -= iv.index.blocksize(is.sector);
        is.sector += 1;
        }
    return is;
    }

QN const& IndexVal::
qn() const 
    { 
    auto is = sectorInfo(*this);
    return index.qn(is.sector);
    }

IndexVal&  IndexVal::
dag() { index.dag(); return *this; }

class IQIndexDat
    {
    public:
    using storage = std::vector<QNInt>;
    using iterator = storage::iterator;
    using const_iterator = storage::const_iterator;
    private:
    storage iq_;
    public:

    IQIndexDat() { }

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

    long
    blocksize(long i) { return iq_[i-1].second; }

    long
    blocksize0(long i) { return iq_[i].second; }

    QN const&
    qn(long i) { return iq_[i-1].first; }

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

#ifdef DEBUG
#define IQINDEX_CHECK_NULL if(pd == 0) throw ITError("IQIndex storage unallocated");
#else
#define IQINDEX_CHECK_NULL
#endif

void
write(std::ostream & s, QNInt const& q)
    {
    write(s,q.first);
    write(s,q.second);
    }

void
read(std::istream & s, QNInt & q)
    {
    read(s,q.first);
    read(s,q.second);
    }

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

long Index::
nblock() const 
    {
    if(not pd) return 0;
    return static_cast<long>(pd->size());
    }

QN const& Index::
qn(long i) const 
    {
    IQINDEX_CHECK_NULL
#ifdef DEBUG
    if(i > nblock())
        {
        Print(nblock());
        Print(i);
        Error("IQIndex::qn arg out of range");
        }
#endif
    return pd->qn(i);
    }

long Index::
blocksize(long i) const 
    {
    IQINDEX_CHECK_NULL
#ifdef DEBUG
    if(i > nblock())
        {
        Print(nblock());
        Print(i);
        Error("Index::blocksize arg out of range");
        }
#endif
    return pd->blocksize(i);
    }

long Index::
blocksize0(long i) const 
    {
    IQINDEX_CHECK_NULL
#ifdef DEBUG
    if(i >= nblock())
        {
        Print(nblock());
        Print(i);
        Error("Index::blocksize0 arg out of range");
        }
#endif
    return pd->blocksize0(i);
    }

void Index::
makeStorage(qnstorage && qi)
    {
    pd = std::make_shared<IQIndexDat>(std::move(qi));
    }

long
totalM(Index::qnstorage const& storage)
    {
    long tm = 0;
    for(auto& iq : storage)
        {
        tm += iq.second;
        }
    return tm;
    }

Index::
Index(qnstorage && ind_qn, 
      TagSet const& ts)
  : Index(totalM(ind_qn),ts)
    { 
    dir_ = Out;
    makeStorage(std::move(ind_qn));
    }

Index::
Index(qnstorage && ind_qn, 
      Arrow dir, 
      TagSet const& ts)
  : Index(totalM(ind_qn),ts)
    { 
    dir_ = dir;
    makeStorage(std::move(ind_qn));
    }

long
QNblock(Index const& I,
        QN const& Q)
    {
    for(auto n : range1(I.nblock()))
        { 
        if(I.qn(n) == Q) return n;
        }
    if(not hasQNs(I)) Error("Index does not contain any QN blocks");
    println("I = ",I);
    println("Q = ",Q);
    Error("Index does not contain given QN block.");
    return 0l;
    }


long
QNblockSize(Index const& I, 
            QN const& Q)
    { 
    return I.blocksize(QNblock(I,Q));
    }

void Index::
write(std::ostream& s) const 
    { 
    if(!bool(*this)) Error("Index::write: Index is default initialized");
    itensor::write(s,primelevel_);
    itensor::write(s,tags_);
    itensor::write(s,id_);
    itensor::write(s,m_);
    itensor::write(s,dir_);
    if(pd) itensor::write(s,*pd);
    else itensor::write(s,IQIndexDat());
    }

Index& Index::
read(std::istream& s)
    {
    itensor::read(s,primelevel_);
    itensor::read(s,tags_);
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
    itensor::read(s,dir_);
    IQIndexDat dat;
    itensor::read(s,dat);
    if(dat.size()>0)
        pd = std::make_shared<IQIndexDat>(std::move(dat.store()));

#ifdef DEBUG
    if(primelevel_ < 0) Error("Negative primeLevel");
#endif

    return *this;
    }

} //namespace itensor

