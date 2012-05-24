//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"
#include <set>
using namespace std;
using boost::format;


void 
intrusive_ptr_add_ref(IQTDat* p) 
    { 
    ++(p->numref); 
    }

void 
intrusive_ptr_release(IQTDat* p) 
    { 
    if(--(p->numref) == 0)
        { 
        delete p; 
        } 
    }

IQTDat::
IQTDat() 
    : 
    numref(0), 
    rmap_init(false)
    { }

IQTDat::
IQTDat(const IQTDat& other) 
    : 
    itensor(other.itensor), 
    numref(0), 
    rmap_init(false)
	{ }

IQTDat::
IQTDat(istream& s) 
    : 
    numref(0),
    rmap_init(false)
    { read(s); }

void IQTDat::
read(istream& s)
	{
	uninit_rmap();
	size_t size;
	s.read((char*) &size,sizeof(size));
	itensor.resize(size);
    Foreach(ITensor& t, itensor)
        { t.read(s); }
	}

void IQTDat::
write(ostream& s) const
	{
	size_t size = itensor.size();
	s.write((char*) &size,sizeof(size));
    Foreach(const ITensor& t, itensor)
        { t.write(s); }
	}

void IQTDat::
init_rmap() const
	{
	if(rmap_init) return;

    for(iterator it = itensor.begin(); it != itensor.end(); ++it)
	    rmap[ApproxReal(it->uniqueReal())] = it;

	rmap_init = true;
	}

void IQTDat::
uninit_rmap() const 
	{ 
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("uninit_rmap called on shared IQTDat");
        }
#endif
	rmap.clear();
	rmap_init = false; 
	}

bool IQTDat::
has_itensor(const ApproxReal& r) const
	{ 
	init_rmap();
	return rmap.count(r) == 1; 
	}

void IQTDat::
clear()
    {
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("clear called on shared IQTDat");
        }
#endif
    uninit_rmap();
    itensor.clear();
    }

void IQTDat::
insert(const ApproxReal& r, const ITensor& t)
    {
    init_rmap();
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("insert called on shared IQTDat");
        }
#endif
    if(rmap.count(r) == 1)
        {
        Print((*rmap[r])); 
        Print(t);
        Error("Can't insert ITensor with identical structure twice, use operator+=.");
        }
    else
        {
        itensor.push_front(t);
        rmap[r] = itensor.begin();
        }
    }

void IQTDat::
insert(const ITensor& t)
    {
    ApproxReal r(t.uniqueReal());
    insert(r,t);
    }

void IQTDat::
insert_add(const ApproxReal& r, const ITensor& t)
    {
    init_rmap();
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("insert_add called on shared IQTDat");
        }
#endif
    if(rmap.count(r) == 1)
        {
        *rmap[r] += t;
        return;
        }
    else
        {
        itensor.push_front(t);
        rmap[r] = itensor.begin();
        }
    }

void IQTDat::
insert_add(const ITensor& t)
    {
    ApproxReal r(t.uniqueReal());
    insert_add(r,t);
    }

void IQTDat::
insert_assign(const ITensor& t)
    {
    init_rmap();
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("insert_assign called on shared IQTDat");
        }
#endif
    ApproxReal r(t.uniqueReal());
    if(rmap.count(r) == 1)
        {
        rmap[r]->assignFrom(t);
        return;
        }
    else
        {
        itensor.push_front(t);
        rmap[r] = itensor.begin();
        }
    }

void IQTDat::
clean(Real min_norm)
    {
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("clean called on shared IQTDat");
        }
#endif
    list<ITensor> nitensor;
    Foreach(const ITensor& t, itensor)
        {
        if(t.norm() >= min_norm)
            nitensor.push_back(t);
        }
    swap(nitensor);
    }

void IQTDat::
swap(StorageT& new_itensor)
    {
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("clean called on shared IQTDat");
        }
#endif
    uninit_rmap();
    itensor.swap(new_itensor);
    }

void IQTDat::
scaleTo(const LogNumber& newscale)
    {
    Foreach(ITensor& t, itensor)
        t.scaleTo(newscale);
    }

//
// IQTensor
//

int IQTensor::
r() const 
    { 
    if(is_ == 0) return 0;
    return is_->r(); 
    }

const IQIndex& IQTensor::
index(int j) const 
    { 
    if(is_ == 0) Error("IQTensor is null");
    return is_->index(j); 
    }

int IQTensor::
iten_size() const { return dat().size(); }

bool IQTensor::
iten_empty() const { return dat().empty(); }

//----------------------------------------------------
//IQTensor: iterators 
IQTensor::const_iten_it IQTensor::
const_iten_begin() const 
    { 
    return dat().begin();
    }

IQTensor::const_iten_it IQTensor::
const_iten_end() const 
    { 
    return dat().end();
    }

std::pair<IQTensor::const_iten_it,IQTensor::const_iten_it> IQTensor::
itensors() const 
    { return std::make_pair(dat().begin(),dat().end()); }

IQIndexSet::index_it IQTensor::
const_iqind_begin() const { return is_->begin(); }

IQIndexSet::index_it IQTensor::
const_iqind_end()   const { return is_->end(); }

std::pair<IQIndexSet::index_it,IQIndexSet::index_it> IQTensor::
iqinds() const 
    { return is_->index(); }

//----------------------------------------------------
//IQTensor: Constructors 

IQTensor::
IQTensor() 
    : 
    is_(0),
    p(0) 
    { }

IQTensor::
IQTensor(Real val) 
    : 
    is_(new IQIndexSet()),
    p(new IQTDat())
    { 
    operator+=(ITensor(val));
    }

IQTensor::
IQTensor(const IQIndex& i1) 
    : 
    is_(new IQIndexSet(i1)),
    p(new IQTDat()) 
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2) 
    : 
    is_(new IQIndexSet(i1,i2)),
    p(new IQTDat()) 
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) 
	: 
    is_(new IQIndexSet(i1,i2,i3)),
    p(new IQTDat()) 
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4) 
    : 
    is_(new IQIndexSet(i1,i2,i3,i4)),
    p(new IQTDat())
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5) 
    : 
    is_(new IQIndexSet(i1,i2,i3,i4,i5)),
    p(new IQTDat()) 
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6)
	: 
    is_(new IQIndexSet(i1,i2,i3,i4,i5,i6)),
    p(new IQTDat()) 
	{ 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7)
	: 
    is_(new IQIndexSet(i1,i2,i3,i4,i5,i6,i7)),
    p(new IQTDat()) 
	{ 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7,const IQIndex& i8)
	: 
    is_(new IQIndexSet(i1,i2,i3,i4,i5,i6,i7,i8)),
    p(new IQTDat()) 
	{ 
    }

IQTensor::
IQTensor(std::vector<IQIndex>& iqinds_) 
	: 
    is_(new IQIndexSet(iqinds_)),
    p(new IQTDat()) 
	{ 
#ifdef DEBUG
    Foreach(const IQIndex& I, iqinds_)
        {
        if(I == IQIndex::Null())
            Error("IQIndex is null");
        }
#endif
    }

IQTensor::
IQTensor(const IQIndexVal& iv1) 
    : 
    is_(new IQIndexSet(iv1.iqind)),
    p(new IQTDat()) 
	{ 
	operator()(iv1) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2) 
	: 
    is_(new IQIndexSet(iv1.iqind,iv2.iqind)),
    p(new IQTDat()) 
	{ 
    operator()(iv1,iv2) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
         const IQIndexVal& iv3) 
	: 
    is_(new IQIndexSet(iv1.iqind,iv2.iqind,iv3.iqind)),
    p(new IQTDat()) 
	{ 
    operator()(iv1,iv2,iv3) = 1;
	}

IQTensor::
IQTensor(ITensor::ITmaker itm) 
    : 
    is_(new IQIndexSet(IQIndex::IndReIm())),
    p(new IQTDat()) 
    {
    if(itm == ITensor::makeComplex_1)      
        {
        operator+=(ITensor::Complex_1());
        }
    else if(itm == ITensor::makeComplex_i) 
        {
        operator+=(ITensor::Complex_i());
        }
    }

IQTensor::
IQTensor(IQmaker i) 
    : 
    is_(new IQIndexSet()),
    p(new IQTDat())
    {
    Index s("sing");
    IQIndex single("single",s,QN());
    is_->addindex(single);
    ITensor st(s,1);
    operator+=(st);
    }

IQTensor::
IQTensor(PrimeType pt,const IQTensor& other) 
    : 
    is_(other.is_),
    p(other.p)
    { doprime(pt); }

IQTensor::
IQTensor(std::istream& s)
    : 
    is_(0),
    p(0) 
    { 
    read(s); 
    }

void IQTensor::
read(std::istream& s)
    {
    bool null_;
    s.read((char*) &null_,sizeof(null_));
    if(null_) 
        { *this = IQTensor(); return; }
    is_ = new IQIndexSet(s);
    p = new IQTDat(s);
    }

void IQTensor::
write(std::ostream& s) const
	{
	bool null_ = isNull();
	s.write((char*) &null_,sizeof(null_));
	if(null_) return;
    is_->write(s);
	dat().write(s);
	}

IQTensor& IQTensor::
operator*=(Real fac) 
    { 
    soloDat();

    if(fac == 0) 
        { 
        ncdat().clear(); 
        return *this; 
        }

    Foreach(ITensor& t, ncdat())
        {
        t *= fac;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator/=(Real fac) 
    { 
    soloDat();

    if(fac == 0) 
        { 
        ncdat().clear(); 
        return *this; 
        }

    Foreach(ITensor& t, ncdat())
        {
        t /= fac;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator*=(const LogNumber& lgnum) 
    { 
    soloDat();

    Foreach(ITensor& t, ncdat())
        {
        t *= lgnum;
        }

    return *this; 
    }

void IQTensor::
insert(const ITensor& t) 
    { 
    if(t.scale().sign() != 0)
        {
        soloDat();
        ncdat().insert(t);
        }
    }

IQTensor& IQTensor::
operator+=(const ITensor& t) 
    { 
    if(t.scale().sign() != 0)
        {
        soloDat();
        ncdat().insert_add(t);
        }
    return *this;
    }

//Non-const element access
Real& IQTensor::
operator()(const IQIndexVal& iv1, const IQIndexVal& iv2,
           const IQIndexVal& iv3, const IQIndexVal& iv4, 
           const IQIndexVal& iv5, const IQIndexVal& iv6,
           const IQIndexVal& iv7, const IQIndexVal& iv8)
	{
    soloDat();
    boost::array<IQIndexVal,NMAX+1> iv 
        = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};

    Real ur = 0; 
    int nn = 0; 
    while(GET(iv,nn+1).iqind != IQIndexVal::Null().iqind) 
        ur += GET(iv,++nn).index().uniqueReal(); 
    if(nn != r()) 
        Error("Wrong number of IQIndexVals provided");
    ApproxReal r(ur);

    if(!dat().has_itensor(r))
        {
        std::vector<Index> indices; 
        indices.reserve(nn);
        for(int j = 1; j <= nn; ++j) 
            {
            if(!hasindex(iv[j].iqind)) 
                Error("IQTensor::operator(): IQIndex not found.");
            indices.push_back(iv[j].index());
            }
        ITensor t(indices);
        ncdat().insert_add(r,t);
        }

    return (ncdat().get(r)).operator()(iv1.blockIndexVal(),
                                       iv2.blockIndexVal(),
                                       iv3.blockIndexVal(),
                                       iv4.blockIndexVal(),
                                       iv5.blockIndexVal(),
                                       iv6.blockIndexVal(),
                                       iv7.blockIndexVal(),
                                       iv8.blockIndexVal());
	}

//const element access
Real IQTensor::
operator()(const IQIndexVal& iv1, const IQIndexVal& iv2,
           const IQIndexVal& iv3, const IQIndexVal& iv4, 
           const IQIndexVal& iv5, const IQIndexVal& iv6,
           const IQIndexVal& iv7, const IQIndexVal& iv8) const
	{
    boost::array<IQIndexVal,NMAX+1> iv 
        = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};

    Real ur = 0; 
    int nn = 0; 
    while(GET(iv,nn+1).iqind != IQIndexVal::Null().iqind) 
        ur += GET(iv,++nn).index().uniqueReal(); 
    if(nn != r()) 
        Error("Wrong number of IQIndexVals provided");
    ApproxReal r(ur);

    if(!dat().has_itensor(r))
        {
        return 0.;
        }
    else
        {
        return (dat().get(r)).operator()(iv1.blockIndexVal(),
                                         iv2.blockIndexVal(),
                                         iv3.blockIndexVal(),
                                         iv4.blockIndexVal(),
                                         iv5.blockIndexVal(),
                                         iv6.blockIndexVal(),
                                         iv7.blockIndexVal(),
                                         iv8.blockIndexVal());
        }
	}

//Method for specifically requesting const access
Real IQTensor::
at(const IQIndexVal& iv1, const IQIndexVal& iv2,
   const IQIndexVal& iv3, const IQIndexVal& iv4, 
   const IQIndexVal& iv5, const IQIndexVal& iv6,
   const IQIndexVal& iv7, const IQIndexVal& iv8) const
    {
    return operator()(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8);
    }


QN IQTensor::
div() const
	{
	QN div_;
	assert(p != 0);
	if(dat().empty())
	    {   
	    this->printIndices("this");
	    Error("IQTensor has no blocks");
	    }
	const ITensor& t = *(dat().begin());
	for(int j = 1; j <= t.r(); ++j)
	    div_ += qn(t.index(j))*dir(t.index(j));
	return div_;
	}


QN IQTensor::
qn(const Index& in) const
	{
	int iqq = find_iqind(in);
	if(iqq == 0) 
	    Error("qn: cant find index");
	return is_->index(iqq).qn(in);
	} 

Arrow IQTensor::
dir(const Index& in) const
	{
	int iqq = find_iqind(in);
	if(iqq == 0) 
	    {
	    this->print("this IQTensor");
	    in.print("in"); 
	    Error("IQTensor::dir(Index&): cant find Index in IQIndices");
	    }
	return is_->index(iqq).dir();
	}

void IQTensor::
ind_inc_prime(const IQIndex& i,int inc)
	{
	soloIndex();
    is_->indIncPrime(i,inc);

    soloDat();

    Foreach(ITensor& t, ncdat())
        {
	    for(int ii = 1; ii <= t.r(); ++ii)
		if(i.hasindex_noprime(t.index(ii)))
		    {
		    int p = t.index(ii).primeLevel();
		    t.mapprimeind(t.index(ii),p,p+inc);
		    }
        }
	}

void IQTensor::
noprime(PrimeType pt)
	{
	solo();

    is_->noprime(pt);

    Foreach(ITensor& t, ncdat())
        { t.noprime(pt); }
	} 

void IQTensor::
noprimelink()
	{
	solo();

    is_->noprime(primeLink);

    Foreach(ITensor& t, ncdat())
        { t.noprime(primeLink); }
	}

void IQTensor::
doprime(PrimeType pt, int inc)
	{
	solo();

    is_->doprime(pt,inc);

    Foreach(ITensor& t, ncdat())
	    t.doprime(pt,inc);
	}

void IQTensor::
mapprime(int plevold, int plevnew, PrimeType pt)
    {
    solo();

    is_->mapprime(plevold,plevnew,pt);

    Foreach(ITensor& t, ncdat())
	    t.mapprime(plevold,plevnew,pt);
	}

void IQTensor::
primeind(const IQIndex& I, int inc)
	{
	solo();

    is_->primeind(I,inc);

    Foreach(ITensor& t, ncdat())
    for(std::vector<inqn>::const_iterator
        x = I.iq().begin(); x != I.iq().end(); ++x)
        {
		if(t.hasindex(x->index)) 
		    t.primeind(x->index,inc);
        }
	}

void IQTensor::
noprimeind(const IQIndex& I)
	{
	solo();

    is_->noprimeind(I);

    Foreach(ITensor& t, ncdat())
    for(std::vector<inqn>::const_iterator
        x = I.iq().begin(); x != I.iq().end(); ++x)
        {
        if(t.hasindex(x->index)) 
            t.noprimeind(x->index);
        }
	}


int IQTensor::
find_iqind(const Index& ii) const
    {
    for(int j = 1; j <= is_->r(); ++j)
        {
        if(is_->index(j).hasindex(ii)) 
            return j;
        }
    return 0;
    }

bool IQTensor::
uses_ind(const Index& ii) const
    {
    Foreach(const ITensor& t, dat())
        {
        if(t.hasindex(ii)) 
            return true;
        }
    return false;
    }

int IQTensor::
findindex(const IQIndex& I) const 
    { 
    if(is_ == 0) return 0;
    return is_->findindex(I); 
    }

bool IQTensor::
hastype(IndexType t) const 
    { 
    if(is_ == 0) return false;
    return is_->hastype(t); 
    }

const IQIndex& IQTensor::
findtype(IndexType t) const 
    { 
    if(is_ == 0) Error("findtype failed, IQTensor is null");
    return is_->findtype(t); 
    }

const IQIndex& IQTensor::
finddir(Arrow dir) const 
    { 
    if(is_ == 0) Error("finddir failed, IQTensor is null");
    return is_->finddir(dir); 
    }

bool IQTensor::
hasindex(const IQIndex& I) const 
    { 
    if(is_ == 0) return false;
    return is_->hasindex(I); 
    }

void IQTensor::
symmetricDiag11(const IQIndex& i1, IQTensor& D, IQTensor& U, IQIndex& mid) const
    {
    int mink,maxk;
    symmetricDiag11(i1,D,U,mid,mink,maxk);
    }

void IQTensor::
symmetricDiag11(const IQIndex& i1, IQTensor& D, IQTensor& U, IQIndex& mid, int& mink, int& maxk) const
    {
    assert(hasindex(i1));
    assert(hasindex(primed(i1)));
    if(r() != 2) Error("symDiag11: rank must be 2");
    //This method is only intended for IQTensors with zero divergence
    DO_IF_DEBUG(checkDiv(*this,QN());)

    vector<inqn> iq(iten_size()); 
    vector<ITensor> UU(iten_size()); 
    vector<ITensor> d(iten_size()); 

    Real minv = std::numeric_limits<Real>::max(),
         maxv = -minv;
    int w = 0;
    int totk = 0;
    int itMink = -1, itMaxk = -1;
    for(IQTensor::const_iten_it it = const_iten_begin(); it != const_iten_end(); ++it)
        {
        Index mi;
        Index act = it->index(1);
        act.noprime();
        it->symmetricDiag11(act,d[w],UU[w],mi,itMink,itMaxk);
        iq[w] = inqn(mi,i1.qn(act));

        if(d[w](mi(itMink)) < minv)
            {
            mink = totk+itMink;
            minv = d[w](mi(itMink));
            }
        if(d[w](mi(itMaxk)) > maxv)
            {
            maxv = d[w](mi(itMaxk));
            maxk = totk+itMaxk;
            }
        totk += mi.m();
        ++w;
        }

    mid = IQIndex((mid.isNull() ? "mid" : mid.rawname()),iq,i1.dir()*Switch);

    //Insert blocks backwards so that mink and maxk
    //label the appropriate values (because IQTDat::insert_add
    //pushes from the front)
    U = IQTensor(i1,mid);
    for(int j = UU.size()-1; j >= 0; --j)
        U.insert(UU.at(j));

    D = IQTensor(mid);
    for(int j = d.size()-1; j >= 0; --j)
        D.insert(d.at(j));
    }

Real IQTensor::
uniqueReal() const 
    { 
    if(is_ == 0) Error("IQTensor is null");
    return is_->uniqueReal(); 
    }

Real IQTensor::
norm() const
    {
    Real res = 0;
    Foreach(const ITensor& t, dat())
        { res += sqr(t.norm()); }

    //Even if the ITensor norms aren't separately too
    //large for Real, their sum may be
    if(res > std::numeric_limits<Real>::max())
        {
        throw TooBigForReal("Norm too large for real in IQTensor::norm()");
        }

    return sqrt(res);
    }

Real IQTensor::
sumels() const
    {
    Real res = 0;
    Foreach(const ITensor& t, dat())
        { res += t.sumels(); }
    return res;
    }

void IQTensor::
scaleOutNorm() const
    {
    Real f = norm();
    Foreach(const ITensor& t, dat())
        t.scaleTo(f);
    }

void IQTensor::
scaleTo(LogNumber newscale) const
    {
    Foreach(const ITensor& t, dat())
        t.scaleTo(newscale);
    }

void IQTensor::
clean(Real min_norm)
    { 
    soloDat(); 
    ncdat().clean(min_norm); 
    }

void IQTensor::
addindex1(const IQIndex& I)
	{
	if(I.m() != 1) 
	    Error("IQTensor::operator*=(IQIndex): IQIndex must have m == 1.");    
	solo(); 
    Foreach(ITensor& t, ncdat())
        { t.addindex1(I.index(1)); }
	is_->addindex(I);
	}

void IQTensor::
tieIndices(const boost::array<IQIndex,NMAX+1>& indices, int niqind, 
           const IQIndex& tied)
    {
    if(niqind < 1) Error("No IQIndices to tie");

    const int nindex = indices[1].nindex();

    boost::intrusive_ptr<IQTDat> np = new IQTDat();
    boost::intrusive_ptr<IQIndexSet> nis_ = new IQIndexSet(tied);

    int nmatched = 0;
    for(int k = 1; k <= is_->r(); ++k)
        {
        const IQIndex& K = is_->index(k);
        bool K_is_tied = false;
        for(int j = 1; j <= niqind; ++j)
        if(K == indices[j]) 
            {
            if(indices[j].m() != tied.m())
                Error("Tied indices must have matching m's");
            K_is_tied = true;
            ++nmatched;
            break;
            }
        if(!K_is_tied)
            {
            nis_->addindex(K);
            }
        }

    if(nmatched != niqind)
        {
        PrintIndices((*this));
        cout << "Indices to tie = " << endl;
        for(int j = 1; j <= niqind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find IQIndex to tie");
        }

    boost::array<Index,NMAX+1> totie;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 1; n <= niqind; ++n)
            totie[n] = indices[n].index(i);

        Foreach(const ITensor& t, dat())
            {
            if(!t.hasAllIndex(totie,niqind)) continue;
            ITensor nt(t);
            nt.tieIndices(totie,niqind,tied.index(i));
            np->insert_add(nt);
            }
        }
    p.swap(np);
    is_.swap(nis_);
    }

void IQTensor::
tieIndices(const IQIndex& i1, const IQIndex& i2, const IQIndex& tied)
    {
    boost::array<IQIndex,NMAX+1> inds =
        {{ IQIndex::Null(), i1, i2, 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null() }};

    tieIndices(inds,2,tied);
    }

void IQTensor::
trace(const boost::array<IQIndex,NMAX+1>& indices, int niqind)
    {
    if(niqind < 1) Error("No IQIndices to trace");

    const int nindex = indices[1].nindex();
    const int tm = indices[1].m();

    boost::intrusive_ptr<IQTDat> np = new IQTDat();
    boost::intrusive_ptr<IQIndexSet> nis_ = new IQIndexSet();

    int nmatched = 0;
    for(int k = 1; k <= is_->r(); ++k)
        {
        const IQIndex& K = is_->index(k);
        bool K_traced = false;
        for(int j = 1; j <= niqind; ++j)
        if(K == indices[j]) 
            {
            if(indices[j].m() != tm)
                Error("Traced indices must have matching m's");
            K_traced = true;
            ++nmatched;
            break;
            }
        if(!K_traced)
            {
            nis_->addindex(K);
            }
        }

    if(nmatched != niqind)
        {
        PrintIndices((*this));
        cout << "Indices to trace = " << endl;
        for(int j = 1; j <= niqind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find IQIndex to trace");
        }

    boost::array<Index,NMAX+1> totrace;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 1; n <= niqind; ++n)
            totrace[n] = indices[n].index(i);

        Foreach(const ITensor& t, dat())
            {
            if(!t.hasAllIndex(totrace,niqind)) continue;
            ITensor tt(t);
            tt.trace(totrace,niqind);
            np->insert_add(tt);
            }
        }
    p.swap(np);
    is_.swap(nis_);
    }

void IQTensor::
trace(const IQIndex& i1, const IQIndex& i2)
    {
    boost::array<IQIndex,NMAX+1> inds =
        {{ IQIndex::Null(), i1, i2, 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null() }};

    trace(inds,2);
    }

void IQTensor::
trace(const IQIndex& i1)
    {
    boost::array<IQIndex,NMAX+1> inds =
        {{ IQIndex::Null(), i1, IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null() }};

    trace(inds,1);
    }

Real
trace(IQTensor T)
    {
    boost::array<IQIndex,NMAX+1> inds;
    for(int k = 1; k <= T.r(); ++k)
        {
        inds.at(k) = T.index(k);
        }
    T.trace(inds,T.r());
    T *= IQTensor::Sing();
    return ReSingVal(T);
    }

int IQTensor::
vecSize() const
	{
    if(this->isNull()) return 0;
	int s = 0;
    Foreach(const ITensor& t, (dat()))
	    s += t.vecSize();
	return s;
	}

int IQTensor::
maxSize() const
	{
    int ms = 1;
	for(int j = 0; j < is_->r(); ++j)
	    ms *= is_->m(j);
    return ms;
    }

void IQTensor::
assignToVec(VectorRef v) const
	{
	if(vecSize() != v.Length())
	    Error("Mismatched sizes in IQTensor::assignToVec(VectorRef v).");
	int off = 1;
	for(const_iten_it jj = const_iten_begin(); jj != const_iten_end(); ++jj)
	    {
	    int d = jj->vecSize();
	    jj->assignToVec(v.SubVector(off,off+d-1));
	    off += d;
	    }
	}

void IQTensor::
assignFromVec(VectorRef v)
	{
	soloDat();
	if(vecSize() != v.Length())
	    Error("bad size");
	int off = 1;
    Foreach(ITensor& t, ncdat())
	    {
	    int d = t.vecSize();
	    t.assignFromVec(v.SubVector(off,off+d-1));
	    off += d;
	    }
	}

void IQTensor::
Randomize() 
	{ 
	soloDat(); 
    Foreach(ITensor& t, ncdat())
        { t.Randomize(); }
	}

void IQTensor::
print(std::string name,Printdat pdat) const 
	{ 
    bool savep = Globals::printdat();
    Globals::printdat() = (pdat==ShowData); 
	std::cerr << "\n" << name << " =\n" << *this << "\n"; 
    Globals::printdat() = savep;
	}

void IQTensor::
printIndices(const std::string& name) const
	{ 
	cout << "\n" << name << " (IQIndices only) = \n";
    if(this->isNull())
        {
        cout << "    [IQTensor is null]" << endl;
        return;
        }
	for(int j = 1; j <= is_->r(); ++j)
	    cout << is_->index(j) << "\n\n";
	cout << "---------------------------\n" << endl;
	}

void IQTensor::
assignFrom(const IQTensor& other)
	{
    //TODO: account for fermion sign here
    if(fabs(uniqueReal()-other.uniqueReal()) > 1E-10)
        {
        PrintIndices((*this));
        PrintIndices(other);
        Error("Can't assign from IQTensor, different IQIndices");
        }

    soloDat();
    Foreach(const ITensor& t, other.dat())
	    {
        ncdat().insert_assign(t);
	    }
	}

void IQTensor::
conj()
    {
    if(!isComplex())
        {
        soloIndex();
        is_->conj();
        return;
        }
    else
        {
        soloDat();

        IQTensor r,i;
        SplitReIm(r,i);

        r.is_->conj();
        i.is_->conj();

        i *= -1.0;
        *this = r * IQTensor::Complex_1() + IQTensor::Complex_i() * i;
        }
    }

void IQTensor::
conj(const IQIndex& I)
    {
    soloIndex();
    is_->conj(I);
    }

void IQTensor::
swap(IQTensor& other)
    {
    is_.swap(other.is_);
    p.swap(other.p);
    }

std::ostream& 
operator<<(std::ostream & s, const IQTensor& T)
    {
    s << "\n----- IQTensor -----\n";
    if(T.isNull())
        {
        s << "(IQTensor is null)\n\n";
        return s;
        }
    s << "IQIndices:\n";
    for(int k = 1; k <= T.r(); ++k)
        { s << "  " << T.index(k) << std::endl; }
    s << "ITensors:\n";
    Foreach(const ITensor& t, T.dat())
        { s << "  " << t << std::endl; }
    s << "-------------------" << "\n\n";
    return s;
    }

void IQTensor::
SplitReIm(IQTensor& re, IQTensor& im) const
    {
    if(!hasindex(IQIndex::IndReIm()))
        {
        IQTensor cop(*this);
        re = cop;
        im = cop;
        im *= 0.0;
        return;
        }
    vector<IQIndex> newreinds;
    remove_copy_if(is_->begin(),is_->end(),std::back_inserter(newreinds),
		    bind2nd(std::equal_to<IQIndex>(),IQIndex::IndReIm()));
    re = IQTensor(newreinds);
    im = re;
    ITensor a,b;
    Foreach(const ITensor& t, dat())
        {
        t.SplitReIm(a,b);
        re.insert(a);
        im.insert(b);
        }
    }

IQTensor& IQTensor::
operator*=(const IQTensor& other)
    {
    //TODO: account for fermion sign here
    if(this == &other)
        {
        IQTensor cp_oth(other);
        return operator*=(cp_oth);
        }

    if(this->isNull()) 
        Error("'This' IQTensor null in product");

    if(other.isNull()) 
        Error("Multiplying by null IQTensor");

    if(hasindex(IQIndex::IndReIm()) && other.hasindex(IQIndex::IndReIm()) && !other.hasindex(IQIndex::IndReImP())
	    && !other.hasindex(IQIndex::IndReImPP()) && !hasindex(IQIndex::IndReImP()) && !hasindex(IQIndex::IndReImPP()))
        {
        static ITensor primer(Index::IndReIm(),Index::IndReImP(),1.0);
        static ITensor primerP(Index::IndReIm(),Index::IndReImPP(),1.0);
        static ITensor prod(Index::IndReIm(),Index::IndReImP(),Index::IndReImPP());
        static IQTensor iqprimer(IQIndex::IndReIm(),IQIndex::IndReImP());
        static IQTensor iqprimerP(IQIndex::IndReIm(),IQIndex::IndReImPP());
        static IQTensor iqprod(IQIndex::IndReIm(),IQIndex::IndReImP(),IQIndex::IndReImPP());
        static int first = 1;
        if(first)
            {
            IndexVal iv0(Index::IndReIm(),1), iv1(Index::IndReImP(),1), iv2(Index::IndReImPP(),1);
            iv0.i = 1; iv1.i = 1; iv2.i = 1; prod(iv0,iv1,iv2) = 1.0;
            iv0.i = 1; iv1.i = 2; iv2.i = 2; prod(iv0,iv1,iv2) = -1.0;
            iv0.i = 2; iv1.i = 2; iv2.i = 1; prod(iv0,iv1,iv2) = 1.0;
            iv0.i = 2; iv1.i = 1; iv2.i = 2; prod(iv0,iv1,iv2) = 1.0;
            first = 0;
            iqprimer += primer;
            iqprimerP += primerP;
            iqprod += prod;
            }
        return *this = (*this * iqprimer) * iqprod * (other * iqprimerP);
        }

    solo();

    set<ApproxReal> common_inds;
    
    //Load iqindex_ with those IQIndex's *not* common to *this and other
    static vector<IQIndex> riqind_holder;
    riqind_holder.resize(0);

    for(int i = 1; i <= is_->r(); ++i)
        {
        const IQIndex& I = is_->index(i);
        const_iqind_it f = find(other.is_->begin(),other.is_->end(),I);
        if(f != other.is_->end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() == I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    this->printIndices("*this");
                    other.printIndices("other");
                    cout << "IQIndex from *this = " << I << endl;
                    cout << "IQIndex from other = " << *f << endl;
                    cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                    throw ArrowError("Incompatible arrow directions in IQTensor::operator*=.");
                    }
            for(size_t n = 0; n < I.iq().size(); ++n) 
                { common_inds.insert(ApproxReal(I.iq()[n].index.uniqueReal())); }

            common_inds.insert(ApproxReal(I.uniqueReal()));
            }
        else 
            { 
            riqind_holder.push_back(I); 
            }
        }

    for(int i = 1; i <= other.is_->r(); ++i)
        {
        const IQIndex& I = other.is_->index(i);
        if(!common_inds.count(ApproxReal(I.uniqueReal())))
            { 
            riqind_holder.push_back(I); 
            }
        }

    is_->swapInds(riqind_holder);

    set<ApproxReal> keys;

    list<ITensor> old_itensor; 
    ncdat().swap(old_itensor);

    //com_this maps the uniqueReal of a set of Index's to be contracted over together
    //to those ITensors in *this.itensor having all Index's in that set
    multimap<ApproxReal,const_iten_it> com_this;
    for(const_iten_it tt = old_itensor.begin(); tt != old_itensor.end(); ++tt)
        {
        Real r = 0.0;
        for(int a = 1; a <= tt->r(); ++a)
            {
            if(common_inds.count(ApproxReal(tt->index(a).uniqueReal())))
                { r += tt->index(a).uniqueReal(); }
            }
        com_this.insert(make_pair(ApproxReal(r),tt));
        keys.insert(ApproxReal(r));
        }

    //com_other is the same as com_this but for other
    multimap<ApproxReal,const_iten_it> com_other;
    for(const_iten_it ot = other.const_iten_begin(); ot != other.const_iten_end(); ++ot)
        {
        Real r = 0.0;
        for(int b = 1; b <= ot->r(); ++b)
            {
            if(common_inds.count(ApproxReal(ot->index(b).uniqueReal())))
                { r += ot->index(b).uniqueReal(); }
            }
        com_other.insert(make_pair(ApproxReal(r),ot));
        keys.insert(ApproxReal(r));
        }

    typedef multimap<ApproxReal,const_iten_it>::iterator mit;
    pair<mit,mit> lrange,rrange;
    ITensor tt;
    for(set<ApproxReal>::iterator k = keys.begin(); k != keys.end(); ++k)
        {
        //Equal range returns the begin and end iterators for the sequence
        //corresponding to multimap[key] as a pair
        lrange = com_this.equal_range(*k);
        rrange = com_other.equal_range(*k);

        //Iterate over all ITensors in *this and other sharing
        //the set of contracted Index's corresponding to k
        for(mit ll = lrange.first; ll != lrange.second; ++ll)
        for(mit rr = rrange.first; rr != rrange.second; ++rr)
            {
            //Multiply the ITensors and add into res
            tt = *(ll->second); tt *= *(rr->second);
            if(tt.scale().sign() != 0)
                ncdat().insert_add(tt);
            }
        }

    return *this;

    } //IQTensor& IQTensor::operator*=(const IQTensor& other)

IQTensor& IQTensor::
operator/=(const IQTensor& other)
    {
    //TODO: account for fermion sign here
    if(this == &other)
        {
        IQTensor cp_oth(other);
        return operator/=(cp_oth);
        }

    if(this->isNull()) 
        Error("'This' IQTensor null in product");

    if(other.isNull()) 
        Error("Multiplying by null IQTensor");

    if(hasindex(IQIndex::IndReIm()) && other.hasindex(IQIndex::IndReIm()) && !other.hasindex(IQIndex::IndReImP())
	    && !other.hasindex(IQIndex::IndReImPP()) && !hasindex(IQIndex::IndReImP()) && !hasindex(IQIndex::IndReImPP()))
        {
        Error("IQTensor::operator/= not yet implemented for complex numbers");
        }


    set<ApproxReal> common_inds;
    
    static vector<IQIndex> riqind_holder;
    riqind_holder.resize(0);

    for(int i = 1; i <= is_->r(); ++i)
        {
        const IQIndex& I = is_->index(i);
        const_iqind_it f = find(other.is_->begin(),other.is_->end(),I);
        if(f != other.is_->end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() != I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    this->printIndices("*this");
                    other.printIndices("other");
                    cout << "IQIndex from *this = " << I << endl;
                    cout << "IQIndex from other = " << *f << endl;
                    cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                    throw ArrowError("Incompatible arrow directions in IQTensor::operator/=.");
                    }
            for(size_t n = 0; n < I.iq().size(); ++n) 
                { common_inds.insert(ApproxReal(I.iq()[n].index.uniqueReal())); }

            common_inds.insert(ApproxReal(I.uniqueReal()));
            }
        riqind_holder.push_back(I);
        }

    bool inds_from_other = false;
    for(int i = 1; i <= other.is_->r(); ++i)
        {
        const IQIndex& I = other.is_->index(i);
        if(!common_inds.count(ApproxReal(I.uniqueReal())))
            { 
            riqind_holder.push_back(I);
            inds_from_other = true;
            }
        }

    //Only update IQIndices if they are different
    //from current set
    if(inds_from_other)
        {
        soloIndex();
        is_->swapInds(riqind_holder);
        }


    soloDat();

    list<ITensor> old_itensor; 
    ncdat().swap(old_itensor);

    set<ApproxReal> keys;

    //com_this maps the uniqueReal of a set of Index's to be summed over together
    //to those ITensors in *this.itensor having all Index's in that set
    multimap<ApproxReal,const_iten_it> com_this;
    for(const_iten_it tt = old_itensor.begin(); tt != old_itensor.end(); ++tt)
        {
        Real r = 0.0;
        for(int a = 1; a <= tt->r(); ++a)
            {
            Real ur = tt->index(a).uniqueReal();
            if(common_inds.count(ApproxReal(ur)))
                r += ur;
            }
        com_this.insert(make_pair(ApproxReal(r),tt));
        keys.insert(ApproxReal(r));
        }

    //com_other is the same as com_this but for other
    multimap<ApproxReal,const_iten_it> com_other;
    for(const_iten_it ot = other.const_iten_begin(); ot != other.const_iten_end(); ++ot)
        {
        Real r = 0.0;
        for(int b = 1; b <= ot->r(); ++b)
            {
            Real ur = ot->index(b).uniqueReal();
            if(common_inds.count(ApproxReal(ur)))
                r += ur;
            }
        com_other.insert(make_pair(ApproxReal(r),ot));
        keys.insert(ApproxReal(r));
        }

    typedef multimap<ApproxReal,const_iten_it>::iterator mit;
    pair<mit,mit> lrange,rrange;
    ITensor tt;
    for(set<ApproxReal>::iterator k = keys.begin(); k != keys.end(); ++k)
        {
        //Equal range returns the begin and end iterators for the sequence
        //corresponding to multimap[key] as a pair
        lrange = com_this.equal_range(*k);
        rrange = com_other.equal_range(*k);

        //Iterate over all ITensors in *this and other sharing
        //the set of contracted Index's corresponding to k
        for(mit ll = lrange.first; ll != lrange.second; ++ll)
        for(mit rr = rrange.first; rr != rrange.second; ++rr)
            {
            //Multiply the ITensors and add into res
            tt = *(ll->second); tt /= *(rr->second);
            if(tt.scale().sign() != 0)
                ncdat().insert_add(tt);
            }
        }

    return *this;

    } //IQTensor& IQTensor::operator/=(const IQTensor& other)

//Extracts the real and imaginary parts of the 
//component of a rank 0 tensor (scalar)
void IQTensor::
GetSingComplex(Real& re, Real& im) const
    {
#ifdef DEBUG
    if(r() != 1)
        {
        PrintIndices((*this));
        Error("GetSingComplex only valid for rank 1 IQTensor with index 'single'");
        }
#endif
    IQTensor tre,tim;
    SplitReIm(tre,tim);

    //Only IQIndex should be IQIndex::IndReIm()
    /*
    if(tre.p->iqindex_.size() != 1)
	{
        cout << *this;
        cout << tre;
        Error("bad tre size");
	}
    if(tim.p->iqindex_.size() != 1) Error("bad tim size");
    */
    for(const_iqind_it jj = tre.is_->begin(); jj != tre.is_->end(); ++jj)
        {
        if(*jj != IQTensor::Sing().is_->index(1))
            {
            cout << *this;
            cout << tre;
            Error("bad tre size");
            }
        }
    for(const_iqind_it jj = tim.is_->begin(); jj != tim.is_->end(); ++jj)
        if(*jj != IQTensor::Sing().is_->index(1))
        { Error("bad tim size"); }

    if(tre.iten_size() == 0)
        { re = 0.0; }
    else
        {
        const ITensor& t = *(tre.dat().begin());
        if(t.vecSize() != 1) 
            {
            cout << "tre is\n" << tre << endl;
            Error("bad tre dat size");
            }
        re = t.val0();
        }
    if(tim.iten_size() == 0)
        { im = 0.0; }
    else
        {
        const ITensor& t = *(tim.dat().begin());
        if(t.vecSize() != 1) Error("bad tim dat size");
        im = t.val0();
        }
    }

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    //TODO: account for fermion sign here
    if(this == &other) 
        {
        operator*=(2);
        return *this;
        }

    IQTensor& This = *this;
    if(fabs(This.uniqueReal()-other.uniqueReal()) > 1.0e-11) 
        {
        PrintIndices(This);
        PrintIndices(other);
        Error("bad match unique real in IQTensor::operator+=");
        }

    if(is_->r() == 0)	// Automatic initializing a summed IQTensor in a loop
        { 
        return (*this = other); 
        }

    bool complex_this = hasindex(IQIndex::IndReIm()); 
    bool complex_other = other.hasindex(IQIndex::IndReIm()); 
    if(!complex_this && complex_other)
        return (This = (This * IQTensor::Complex_1()) + other);
    if(complex_this && !complex_other)
        return (This += other * IQTensor::Complex_1());

    soloDat(); 

    Foreach(const ITensor& t, other.dat())
        { 
        ncdat().insert_add(t);
        }

    return *this;
    }

//
//Automatically convert this IQTensor
//to an ITensor
//
ITensor IQTensor::
toITensor() const
    {
    //Resulting ITensor's indices are 
    //the Index versions of this's IQIndices
    vector<Index> indices;
    for(int j = 1; j <= is_->r(); ++j)
        {
        indices.push_back(Index(is_->index(j)));
        }
    ITensor res(indices);

    //Loop over ITensors (blocks) within this IQTensor
    Foreach(const ITensor& t, dat())
        {
        ITensor exp(t);
        //Loop over Index's of the k'th ITensor
        for(int j = 1; j <= t.r(); ++j)
            {
            //Want to transform 'small' into the 
            //Index version of the IQIndex that contains
            //it, with the appropriate offset
            const Index& small = t.index(j);
            //Find the IQIndex that contains 'small'
            const IQIndex* big = 0;
            int offset = -1;
            for(int q = 1; q <= is_->r(); ++q)
                if(is_->index(q).hasindex(small))
                    {
                    big = &(is_->index(q));
                    offset = big->offset(small);
                    break;
                    }
            exp.expandIndex(small,*big,offset);
            }
        //Once all Indices expanded, add to res
        res += exp;
        }
    return res;
    } //IQTensor::operator ITensor() const

void IQTensor::
soloIndex()
	{
	if(is_ == 0)
        Error("IQTensor is null");

	if(is_->count() != 1)
        is_ = new IQIndexSet(*is_);
    }

void IQTensor::
soloDat()
    {
	if(p == 0)
        Error("IQTensor is null");

	if(dat().count() != 1)
        p = new IQTDat(*p);
	}

void IQTensor::
solo()
    {
    soloIndex();
    soloDat();
    }

Real 
ReSingVal(const IQTensor& x)
    {
    Real re, im;
    x.GetSingComplex(re,im);
    return re;
    }

Real 
Dot(const IQTensor& x, const IQTensor& y)
    {
    IQTensor res(y);

    int j = x.findindex(y.index(1));
    if(x.index(j).dir() != y.index(1).dir())
        {
        //Arrows do match
        res *= x;
        }
    else
        {
        //Arrows do not match, fix:
        res *= conj(x);
        }
    res *= IQTensor::Sing();
    return ReSingVal(res);
    }

void 
BraKet(const IQTensor& x, const IQTensor& y, Real& re, Real& im)
    {
    IQTensor res(IQTensor::Sing() * conj(x) * y);
    res.GetSingComplex(re,im);
    }

void 
checkQNs(const IQTensor& T)
    {
    QN qtot = T.div();
    for(IQTensor::const_iten_it 
        it = T.const_iten_begin(); it != T.const_iten_end(); ++it)
        {
        QN q;
        for(int j = 1; j <= it->r(); ++j) 
            q += T.qn(it->index(j))*T.dir(it->index(j));

        if(q != qtot) 
            {
            std::cout << "\nqtot = " << qtot << "\n\n";
            std::cout << "Offending ITensor = " << *it << "\n\n";
            T.printIndices("T");
            Error("checkQNs: inconsistent QN");
            }
        }
    }

void 
checkDiv(const IQTensor& T, QN expected)
	{
	if(T.dat().empty())
	    {   
	    T.printIndices("this");
	    Error("IQTensor has no blocks");
	    }

    Foreach(const ITensor& t, T.dat())
	    {
	    QN div_;
	    for(int j = 1; j <= t.r(); ++j)
            div_ += T.qn(t.index(j))*T.dir(t.index(j));
	    if(div_ != expected)
            {
            Print(expected);
            Print(div_);
            T.printIndices("This IQTensor:");
            cout << "Incorrect block:\n";
            Print(t);
            Error("Block didn't match expected div");
            }
	    }
	}
