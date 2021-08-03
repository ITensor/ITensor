//
// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
#include "itensor/index.h"
#include "itensor/util/readwrite.h"
#include "itensor/util/print_macro.h"

using std::string;
using std::stringstream;
using std::vector;

namespace itensor {


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
    dim_(1),
    tags_(TagSet("0"))
    {
    if(primeLevel() < 0) setPrime(0);
    }

Index::
Index(long m, 
      TagSet const& t)
  : id_(generateID()),
    dim_(m),
    tags_(t)
    { 
    if(primeLevel() < 0) setPrime(0);
    } 

Index::
Index(id_type id,
      long dim, 
      Arrow dir, 
      TagSet const& ts)
  : id_(id),
    dim_(dim),
    dir_(dir),
    tags_(ts)
    { 
    if(primeLevel() < 0) setPrime(0);
    } 

Index::
Index(id_type id,
      long dim, 
      Arrow dir, 
      TagSet const& ts,
      qnstorage&& qns)
  : id_(id),
    dim_(dim),
    dir_(dir),
    tags_(ts)
    { 
    makeStorage(std::move(qns));
    if(primeLevel() < 0) setPrime(0);
    } 


Index& Index::
setPrime(int plev) 
    { 
    tags_.setPrime(plev);
#ifdef DEBUG
    if(this->primeLevel() < 0)
        Error("Negative primeLevel");
#endif
    return *this;
    }


Index& Index::
noPrime()  
    {
    tags_.setPrime(0);
    return *this;
    }


Index& Index::
prime(int inc) 
    { 
    tags_.prime(inc);
#ifdef DEBUG
    if(this->primeLevel() < 0)
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
    return (id(i1) == id(i2)) && (tags(i1) == tags(i2));
    }

bool 
operator!=(Index const& i1, Index const& i2)
    { 
    return not operator==(i1,i2);
    }

Index::id_type
id(Index const& I) { return I.id(); }

long
dim(Index const& I) { return I.dim(); }

long
dim(IndexVal const& I) { return I.dim(); }

int
primeLevel(Index const& I) { return I.primeLevel(); }

TagSet
tags(const Index& I) { return I.tags(); }

Index
addTags(Index I, TagSet const& t) { I.addTags(t); return I; }

Index
removeTags(Index I, TagSet const& t) { I.removeTags(t); return I; }

Index
setTags(Index I, TagSet const& t) { I.setTags(t); return I; }

Index
noTags(Index I) { I.noTags(); return I; }

Index
replaceTags(Index I, TagSet const& tsold, TagSet const& tsnew) { I.replaceTags(tsold,tsnew); return I; }

Index
prime(Index I, int plinc) { I.prime(plinc); return I; }

Index
setPrime(Index I, int plev) { I.setPrime(plev); return I; }

Index
noPrime(Index I) { I.noPrime(); return I; }

IndexVal
addTags(IndexVal I, TagSet const& t) { I.addTags(t); return I; }

IndexVal
removeTags(IndexVal I, TagSet const& t) { I.removeTags(t); return I; }

IndexVal
setTags(IndexVal I, TagSet const& t) { I.setTags(t); return I; }

IndexVal
noTags(IndexVal I) { I.noTags(); return I; }

IndexVal
replaceTags(IndexVal I, TagSet const& tsold, TagSet const& tsnew) { I.replaceTags(tsold,tsnew); return I; }

IndexVal
prime(IndexVal I, int plinc) { I.prime(plinc); return I; }

IndexVal
setPrime(IndexVal I, int plev) { I.setPrime(plev); return I; }

IndexVal
noPrime(IndexVal I) { I.noPrime(); return I; }

//
// Check if Index I contains the tags tsmatch.
//
bool
hasTags(Index I, const TagSet& tsmatch) { return hasTags(tags(I),tsmatch); }

bool
hasQNs(Index const& I) { return nblock(I)!=0; }

Index
removeQNs(Index I) { if(hasQNs(I)) I.removeQNs(); return I; }

Index const&
index(IndexVal const& res) { return res.index; }

long
val(IndexVal const& res) { return res.val; }

bool
hasQNs(IndexVal const& iv) { return hasQNs(index(iv)); }

Index
dag(Index res) { res.dag(); return res; }

IndexVal
dag(IndexVal res) { res.dag(); return res; }

QN const&
qn(IndexVal const& iv) { return iv.qn(); }

Arrow
dir(Index const& res) { return res.dir(); }

Arrow
dir(IndexVal const& res) { return dir(index(res)); }

//TODO: what is this for? It doesn't make as much sense with
//tags, but I guess we can compare tags with inequalities
bool
operator>(Index const& i1, Index const& i2)
    { 
    if(dim(i1) == dim(i2)) 
        {
        if(id(i1) == id(i2)) return primeLevel(i1) > primeLevel(i2);
        return id(i1) > id(i2);
        }
    return dim(i1) > dim(i2);
    }

bool
operator<(Index const& i1, Index const& i2)
    {
    if(dim(i1) == dim(i2)) 
        {
        if(id(i1) == id(i2)) return primeLevel(i1) < primeLevel(i2);
        return id(i1) < id(i2);
        }
    return dim(i1) < dim(i2);
    }

std::ostream& 
operator<<(std::ostream & s, Index const& I)
    {
    s << "(dim=";
    s << dim(I);
    if(Global::showIDs()) 
        {
        s << "|id=" << (id(I) % 1000);
        //s << "," << id(I);
        }
    auto ts = tags(I);
    if(size(ts) > 0)
      {
      s << "|\"";
      auto ts = tags(I);
      for(auto i : range(size(ts)))
        {
        s << ts[i];
        if( i < (size(ts)-1) ) s << ",";
        }
      s << "\"";
      }
    s << ")"; 
    if(primeLevel(I) > 0) 
        {
        if(primeLevel(I) > 3)
            {
            s << "'" << primeLevel(I);
            }
        else
            {
            for(int n = 1; n <= primeLevel(I); ++n)
                s << "'";
            }
        }
    else if(primeLevel(I) < 0) 
        {
        s << "'" << primeLevel(I) << " (WARNING: prime level of Index is negative, not well defined)";
        }
    if(hasQNs(I))
        {
        s << " <" << I.dir() << ">\n";
        for(auto j : range1(nblock(I)))
            {
            s << "  " << j << ": " << blocksize(I,j) << " " <<  qn(I,j);
            if(j != nblock(I)) s << "\n";
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
    //if(val_ < 1 || val_ > dim(index))
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
    return (index(iv1) == index(iv2) && iv1.val == iv2.val);
    }

bool
operator!=(IndexVal const& iv1, IndexVal const& iv2)
    {
    return not operator==(iv1,iv2);
    }

bool
operator==(Index const& I, IndexVal const& iv)
    {
    return index(iv) == I;
    }

bool
operator==(IndexVal const& iv, Index const& I)
    {
    return index(iv) == I;
    }

string
showDim(Index const& I) { return format("dim=%d",dim(I)); }

std::ostream& 
operator<<(std::ostream& s, IndexVal const& iv)
    { 
    return s << "IndexVal: val = " << iv.val 
             << ", ind = " << index(iv);
    }

void
add(Args            & args, 
    Args::Name const& name, 
    TagSet     const& ts) 
    { 
    args.add(name,ts); 
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
    while(is.sind > blocksize(index(iv),is.sector))
        {
        is.sind -= blocksize(index(iv),is.sector);
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

Index
directSum(Index const& i,
          Index const& j,
          Args const& args)
  {
  auto tags = getTagSet(args,"Tags","sum");
  if(not hasQNs(i) && not hasQNs(j))
    {
    auto dim_ij = dim(i) + dim(j);
    if(dim_ij <= 0) dim_ij = 1;
    return Index(dim_ij,tags);
    }
  else
    {
#ifdef DEBUG
    if( dir(i) != dir(j) ) Error("In directSum(Index, Index), input indices must have same arrow direction");
#endif
    auto nblock_i = nblock(i);
    auto nblock_j = nblock(j);
    auto siq = stdx::reserve_vector<QNInt>(nblock_i+nblock_j);
    for(auto iqn : range1(nblock_i))
        siq.emplace_back(qn(i,iqn),blocksize(i,iqn));
    for(auto jqn : range1(nblock_j))
        siq.emplace_back(qn(j,jqn),blocksize(j,jqn));
#ifdef DEBUG
    if(siq.empty()) Error("siq is empty in plussers");
#endif
    return Index(std::move(siq),dir(i),tags);
    }
  return Index();
  }

Index
sim(Index const& I)
    {
    Index J;
    J.sim(I);
    return J;
    }

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

long
nblock(Index const& i) { return i.nblock(); }

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

QN const&
qn(Index const& i, long b) { return i.qn(b); }

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

long
blocksize(Index const& i, long b) { return i.blocksize(b); }

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
    if(primeLevel() < 0) setPrime(0);
    }

Index::
Index(qnstorage && ind_qn, 
      Arrow dir, 
      TagSet const& ts)
  : Index(totalM(ind_qn),ts)
    { 
    dir_ = dir;
    makeStorage(std::move(ind_qn));
    if(primeLevel() < 0) setPrime(0);
    }

long
QNblock(Index const& I,
        QN const& Q)
    {
    for(auto n : range1(nblock(I)))
        { 
        if(qn(I,n) == Q) return n;
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
    return blocksize(I,QNblock(I,Q));
    }

void Index::
write(std::ostream& s) const 
    { 
    if(!bool(*this)) Error("Index::write: Index is default initialized");
    itensor::write(s,tags_);
    itensor::write(s,id_);
    itensor::write(s,dim_);
    itensor::write(s,dir_);
    if(pd) itensor::write(s,*pd);
    else itensor::write(s,IQIndexDat());
    }

Index& Index::
read(std::istream& s)
    {
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
    itensor::read(s,dim_);
    itensor::read(s,dir_);
    IQIndexDat dat;
    itensor::read(s,dat);
    if(dat.size()>0)
        pd = std::make_shared<IQIndexDat>(std::move(dat.store()));

#ifdef DEBUG
    if(tags_.primeLevel() < 0) Error("Negative primeLevel");
#endif

    return *this;
    }

bool
isFermionic(Index const& I)
    {
    for(auto i : range1(I.nblock()))
        {
        for(auto& qnum : I.qn(i).store())
            {
            if(isFermionic(qnum)) return true;
            }
        }
    return false;
    }

#ifdef ITENSOR_USE_HDF5

void
h5_write(h5::group parent, std::string const& name, Index const& I)
    {
    auto g = parent.create_group(name);
    h5_write_attribute(g,"type","Index",true);
    h5_write_attribute(g,"version",long(1));
    h5_write(g,"id",static_cast<unsigned long>(I.id()));
    h5_write(g,"dim",long(I.dim()));
    h5_write(g,"dir",long(I.dir()));
    h5_write(g,"plev",long(I.primeLevel()));
    h5_write(g,"tags",I.tags());
    if(hasQNs(I))
        {
        h5_write_attribute(g,"space_type","QNBlocks");
        //Write "QNBlocks" data of this QN Index
            {
            auto qg = g.create_group("space");
            h5_write_attribute(qg,"type","QNBlocks",true);
            h5_write_attribute(qg,"version",long(1));
            h5_write(qg,"length",I.nblock());
            auto dims = std::vector<long>(I.nblock());
            for(auto n : range1(I.nblock()))
                {
                dims[n-1] = I.blocksize(n);
                h5_write(qg,format("QN[%d]",n),I.qn(n));
                }
            h5_write(qg,"dims",dims);
            }
        }
    else
        {
        h5_write_attribute(g,"space_type","Int");
        }
    }

void
h5_read(h5::group parent, std::string const& name, Index & I)
    {
    auto g = parent.open_group(name);
    auto type = h5_read_attribute<string>(g,"type");
    if(type != "Index") Error("Group does not contain Index data in HDF5 file");
    auto id = h5_read<unsigned long>(g,"id");
    auto dim = h5_read<long>(g,"dim");
    auto dir = h5_read<long>(g,"dir");
    auto tags = h5_read<TagSet>(g,"tags");
    auto space_type = h5_read_attribute<string>(g,"space_type");
    if(space_type == "QNBlocks") // is a QN Index
        {
        auto qg = g.open_group("space");
        auto nblocks = h5_read<long>(qg,"length");
        auto dims = h5_read<vector<long>>(qg,"dims");
        auto qns = vector<QNInt>(nblocks);
        for(auto n : range1(nblocks))
            {
            auto qn = h5_read<QN>(qg,format("QN[%d]",n));
            qns[n-1] = std::make_pair(qn,dims[n-1]);
            }
        I = Index(id,dim,toArrow(dir),tags,std::move(qns));
        }
    else
        {
        I = Index(id,dim,toArrow(dir),tags);
        }
    auto plev = h5_read<long>(g,"plev");
    I.setPrime(plev);
    }

#endif


} //namespace itensor

