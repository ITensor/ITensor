#include "iqtensor.h"
#include <set>
using namespace std;

DatAllocator<IQTDat> IQTDat::allocator;

IQTDat::
IQTDat() : numref(0), rmap_init(false) { }

IQTDat::
IQTDat(const IQIndex& i1) 
    : iqindex_(1), numref(0), rmap_init(false)
	{ iqindex_[0] = i1; }

IQTDat::
IQTDat(const IQIndex& i1, const IQIndex& i2)
    : iqindex_(2), numref(0), rmap_init(false)
	{ 
	iqindex_[0] = i1; 
	iqindex_[1] = i2; 
	}

IQTDat::
IQTDat(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3)
    : iqindex_(3), numref(0), rmap_init(false)
	{ 
	iqindex_[0] = i1; 
	iqindex_[1] = i2; 
	iqindex_[2] = i3; 
	}

IQTDat::
IQTDat(const IQIndex& i1, const IQIndex& i2, 
       const IQIndex& i3, const IQIndex& i4,
       const IQIndex& i5, const IQIndex& i6, 
	   const IQIndex& i7, const IQIndex& i8)
    : iqindex_(4), numref(0), rmap_init(false)
	{ 
	iqindex_[0] = i1; 
	iqindex_[1] = i2; 
	iqindex_[2] = i3; 
	iqindex_[3] = i4; 
	if(i5 != IQIndex::Null()) 
	    iqindex_.push_back(i5);
	if(i6 != IQIndex::Null()) 
	    iqindex_.push_back(i6);
	if(i7 != IQIndex::Null()) 
	    iqindex_.push_back(i7);
	if(i8 != IQIndex::Null()) 
	    iqindex_.push_back(i8);
	}

IQTDat::
IQTDat(vector<IQIndex>& iqinds_) 
    : numref(0), rmap_init(false) 
    { iqindex_.swap(iqinds_); }

IQTDat::
IQTDat(const IQTDat& other) 
    : itensor(other.itensor), iqindex_(other.iqindex_), numref(0), rmap_init(false)
	{ }

IQTDat::
IQTDat(istream& s) 
    : numref(0) 
    { read(s); }

void IQTDat::
read(istream& s)
	{
	uninit_rmap();
	size_t size;
	s.read((char*) &size,sizeof(size));
	itensor.resize(size);
    for(iten_it it = itensor.begin(); it != itensor.end(); ++it)
        { it->read(s); }

	s.read((char*) &size,sizeof(size));
	iqindex_.resize(size);
    for(iqind_it jj = iqindex_.begin(); jj != iqindex_.end(); ++jj)
        { jj->read(s); }
	}

void IQTDat::
write(ostream& s) const
	{
	size_t size = itensor.size();
	s.write((char*) &size,sizeof(size));
    for(const_iten_it it = itensor.begin(); it != itensor.end(); ++it)
        { it->write(s); }

	size = iqindex_.size();
	s.write((char*) &size,sizeof(size));
    for(const_iqind_it jj = iqindex_.begin(); jj != iqindex_.end(); ++jj)
        { jj->write(s); }
	}

void IQTDat::
init_rmap() const
	{
	if(rmap_init) return;
	for(iten_it it = itensor.begin(); 
        it != itensor.end(); 
        ++it)
	    rmap[ApproxReal(it->uniqueReal())] = it;
	rmap_init = true;
	}

void IQTDat::
uninit_rmap() const 
	{ 
	assert(numref <= 1); 
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
insert_itensor(const ApproxReal& r, const ITensor& t)
    {
    itensor.push_front(t);
    rmap[r] = itensor.begin();
    }

void IQTDat::
clean(Real min_norm)
    {
    list<ITensor> nitensor;
    for(const_iten_it it = itensor.begin(); it != itensor.end(); ++it)
        {
        if(it->norm() >= min_norm)
            nitensor.push_back(*it);
        }
    itensor.swap(nitensor);
    }

int IQTensor::
r() const { return p->iqindex_.size(); }

const IQIndex& IQTensor::
index(int j) const { return GET(p->iqindex_,j-1); }

int IQTensor::
iten_size() const { return p->itensor.size(); }

bool IQTensor::
iten_empty() const { return p->itensor.empty(); }

int IQTensor::
num_index() const { return p->iqindex_.size(); }

//----------------------------------------------------
//IQTensor: iterators 
IQTensor::const_iten_it IQTensor::
const_iten_begin() const { return p->itensor.begin(); }

IQTensor::const_iten_it IQTensor::
const_iten_end() const { return p->itensor.end(); }

std::pair<IQTensor::const_iten_it,IQTensor::const_iten_it> IQTensor::
itensors() const 
    { return std::make_pair(p->itensor.begin(),p->itensor.end()); }

IQTensor::const_iqind_it IQTensor::
const_iqind_begin() const { return p->iqindex_.begin(); }

IQTensor::const_iqind_it IQTensor::
const_iqind_end()   const { return p->iqindex_.end(); }

std::pair<IQTensor::const_iqind_it,IQTensor::const_iqind_it> IQTensor::
iqinds() const 
    { return std::make_pair(p->iqindex_.begin(),p->iqindex_.end()); }

//----------------------------------------------------
//IQTensor: Constructors 

IQTensor::
IQTensor() 
    : p(0) 
    { }

IQTensor::
IQTensor(Real val) 
    : p(new IQTDat()) 
    { 
    operator+=(ITensor(val));
    }

IQTensor::
IQTensor(const IQIndex& i1) 
    : p(new IQTDat(i1)) 
    { }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2) 
    : p(new IQTDat(i1,i2)) 
    { }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) 
	: p(new IQTDat(i1,i2,i3)) 
    { }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4) 
    : p(new IQTDat(i1,i2,i3,i4))
    { }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5) 
    : p(new IQTDat(i1,i2,i3,i4,i5))
    { }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6)
	: p(new IQTDat(i1,i2,i3,i4,i5,i6))
	{ }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7)
	: p(new IQTDat(i1,i2,i3,i4,i5,i6,i7))
	{ }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7,const IQIndex& i8)
	: p(new IQTDat(i1,i2,i3,i4,i5,i6,i7,i8))
	{ }

IQTensor::
IQTensor(std::vector<IQIndex>& iqinds_) 
	: p(new IQTDat(iqinds_))
	{ }

IQTensor::
IQTensor(const IQIndexVal& iv1) 
    : p(new IQTDat(iv1.iqind))
	{ 
	operator()(iv1) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2) 
	: p(new IQTDat(iv1.iqind,iv2.iqind))
	{ 
    operator()(iv1,iv2) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
         const IQIndexVal& iv3) 
	: p(new IQTDat(iv1.iqind,iv2.iqind,iv3.iqind))
	{ 
    operator()(iv1,iv2,iv3) = 1;
	}

IQTensor::
IQTensor(ITensor::ITmaker itm) 
    : p(new IQTDat(IQIndex::IndReIm()))
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
    : p(new IQTDat())
    {
    Index s("sing");
    IQIndex single("single",s,QN());
    p->iqindex_.push_back(single);
    ITensor st(s,1);
    operator+=(st);
    }

IQTensor::
IQTensor(PrimeType pt,const IQTensor& other) 
    : p(other.p)
    { doprime(pt); }

IQTensor::
IQTensor(std::istream& s)
    : p(0) 
    { read(s); }

void IQTensor::
read(std::istream& s)
    {
    bool null_;
    s.read((char*) &null_,sizeof(null_));
    if(null_) 
        { *this = IQTensor(); return; }
    p = new IQTDat(s);
    }

void IQTensor::
write(std::ostream& s) const
	{
	bool null_ = isNull();
	s.write((char*) &null_,sizeof(null_));
	if(null_) return;
	p->write(s);
	}

IQTensor& IQTensor::
operator*=(Real fac) 
    { 
    solo();
    if(fac == 0) 
        { p->itensor.clear(); p->uninit_rmap(); return *this; }
    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        {
        (*it) *= fac;
        }
    return *this; 
    }

void IQTensor::
insert(const ITensor& t) 
	{ 
	solo();
	ApproxReal r(t.uniqueReal());
	if(p->has_itensor(r))
	    {
	    Print(*(p->rmap[r])); 
        Print(t);
	    Error("Can't insert ITensor with identical structure twice, use operator+=.");
	    }
	p->insert_itensor(r,t);
	}

IQTensor& IQTensor::
operator+=(const ITensor& t) 
    { 
    solo();
    ApproxReal r(t.uniqueReal());

    if(t.scale().isRealZero()) { return *this; }

    if(!p->has_itensor(r)) 
        {
        p->insert_itensor(r,t);
        }
    else 
        {
        *(p->rmap[r]) += t;
        }
    return *this;
    }

Real& IQTensor::
operator()(const IQIndexVal& iv1, const IQIndexVal& iv2,
           const IQIndexVal& iv3, const IQIndexVal& iv4, 
           const IQIndexVal& iv5, const IQIndexVal& iv6,
           const IQIndexVal& iv7, const IQIndexVal& iv8)
	{
    solo();
    boost::array<IQIndexVal,NMAX+1> iv 
        = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};

    Real ur = 0; 
    int nn = 0; 
    while(GET(iv,nn+1).iqind != IQIndexVal::Null().iqind) 
        ur += GET(iv,++nn).index().uniqueReal(); 
    if(nn != r()) 
        Error("Wrong number of IQIndexVals provided");
    ApproxReal r(ur);

    if(!p->has_itensor(r))
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
        p->insert_itensor(r,t);
        }
    return (p->rmap[r])->operator()(iv1.blockIndexVal(),
                                    iv2.blockIndexVal(),
                                    iv3.blockIndexVal(),
                                    iv4.blockIndexVal(),
                                    iv5.blockIndexVal(),
                                    iv6.blockIndexVal(),
                                    iv7.blockIndexVal(),
                                    iv8.blockIndexVal());
	}

QN IQTensor::
div() const
	{
	QN div_;
	assert(p != 0);
	if(p->itensor.empty())
	    {   
	    this->printIndices("this");
	    Error("IQTensor has no blocks");
	    }
	const ITensor& t = p->itensor.front();
	for(int j = 1; j <= t.r(); ++j)
	    div_ += qn(t.index(j))*dir(t.index(j));
	return div_;
	}

void IQTensor::
checkDiv(QN expected) const
	{
	assert(p != 0);
	if(p->itensor.empty())
	    {   
	    this->printIndices("this");
	    Error("IQTensor has no blocks");
	    }
    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
	    {
        const ITensor& t = *it;
	    QN div_;
	    for(int j = 1; j <= t.r(); ++j)
		div_ += qn(t.index(j))*dir(t.index(j));
	    if(div_ != expected)
            {
            Print(expected);
            Print(div_);
            this->printIndices("this IQTensor:");
            std::cout << "Incorrect block:\n";
            Print(t);
            Error("Block didn't match expected div");
            }
	    }
	}

QN IQTensor::
qn(const Index& in) const
	{
	int iqq = find_iqind(in)-1;
	if(iqq == -1) 
	    Error("qn: cant find index");
	return p->iqindex_[iqq].qn(in);
	} 

Arrow IQTensor::
dir(const Index& in) const
	{
	int iqq = find_iqind(in)-1;
	if(iqq == -1) 
	    {
	    this->print("this IQTensor");
	    in.print("in"); 
	    Error("IQTensor::dir(Index&): cant find Index in IQIndices");
	    }
	return p->iqindex_[iqq].dir();
	}

void IQTensor::
ind_inc_prime(const IQIndex& i,int inc)
	{
	solo();
	p->uninit_rmap();
	bool gotit = false;
    for(iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
	    if(jj->noprime_equals(i))
		{
		gotit = true;
		int p = jj->primeLevel();
		jj->mapprime(p,p+inc);
		}

	if(!gotit)
	    {
	    std::cerr << "IQIndex was " << i << "\n";
	    std::cout << "IQIndex was " << i << "\n";
	    Error("ind_inc_prime: couldn't find IQIndex");
	    }

	for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    for(int ii = 1; ii <= jj->r(); ++ii)
		if(i.hasindex_noprime(jj->index(ii)))
		    {
		    int p = jj->index(ii).primeLevel();
		    jj->mapprimeind(jj->index(ii),p,p+inc);
		    }
	}

void IQTensor::
noprime(PrimeType pt)
	{
	solo();
	p->uninit_rmap();
    for(iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        { jj->noprime(pt); }
	for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        { it->noprime(pt); }
	} 

void IQTensor::
noprimelink()
	{
	solo();
	p->uninit_rmap();
    for(iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        {
	    if(jj->type() == Link) 
		jj->noprime();
        }
	for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        { it->noprime(primeLink); }
	}

void IQTensor::
doprime(PrimeType pt, int inc)
	{
	solo();
	p->uninit_rmap();

	DoPrimer prim(pt,inc);
	for_each(p->iqindex_.begin(), p->iqindex_.end(),prim);

	for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    jj->doprime(pt,inc);
	}

void IQTensor::
mapprime(int plevold, int plevnew, PrimeType pt)
    {
    solo();
    p->uninit_rmap();

	MapPrimer prim(plevold,plevnew,pt);
	for_each(p->iqindex_.begin(), p->iqindex_.end(),prim);

	for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    jj->mapprime(plevold,plevnew,pt);
	}

void IQTensor::
primeind(const IQIndex& I)
	{
	solo();
	p->uninit_rmap();
    for(iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        {
	    if(*jj == I) 
            {
            *jj = primed(*jj);
            break;
            }
        }

    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
    for(std::vector<inqn>::const_iterator
        x = I.iq().begin(); x != I.iq().end(); ++x)
        {
		if(it->hasindex(x->index)) 
		    it->primeind(x->index);
        }
	}

void IQTensor::
noprimeind(const IQIndex& I)
	{
	solo();
	p->uninit_rmap();
    for(iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        {
	    if(*jj == I) 
            {
            jj->noprime();
            break;
            }
        }

    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
    for(std::vector<inqn>::const_iterator
        x = I.iq().begin(); x != I.iq().end(); ++x)
        {
        if(it->hasindex(x->index)) 
            it->noprimeind(x->index);
        }
	}


int IQTensor::
find_iqind(const Index& I) const
    {
    for(size_t j = 0; j < p->iqindex_.size(); ++j)
        {
        if(p->iqindex_[j].hasindex(I)) 
            return j+1;
        }
    return 0;
    }

bool IQTensor::
uses_ind(const Index& i) const
    {
    for(const_iten_it it=p->itensor.begin(); it != p->itensor.end(); ++it)
        {
        if(it->hasindex(i)) 
            return true;
        }
    return false;
    }

int IQTensor::
findindex(const IQIndex& i) const
    {
    const_iqind_it f = find(p->iqindex_.begin(),p->iqindex_.end(),i);
    if(f == p->iqindex_.end()) 
        return 0;
    else 
        return (f - p->iqindex_.begin())+1;
    }

bool IQTensor::
hastype(IndexType t) const
    {
    for(const_iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        { if(jj->type() == t) return true; }
    return false;
    }

const IQIndex& IQTensor::
findtype(IndexType t) const
    {
    for(const_iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        {
        if(jj->type() == t) 
        return *jj;
        }
    Error("IQTensor::findtype: couldn't find type");
    return IQIndex::Null();
    }

const IQIndex& IQTensor::
finddir(Arrow dir) const
    {
    for(const_iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        { if(jj->dir() == dir) return *jj; }
    Error("IQTensor::finddir: couldn't find dir");
    return IQIndex::Null();
    }

bool IQTensor::
hasindex(const IQIndex& i) const 
    { 
    for(size_t j = 0; j < p->iqindex_.size(); ++j)
        if(i == p->iqindex_[j]) 
        return true;
    return false;
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
    DO_IF_DEBUG(this->checkDiv(QN());)

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
    //label the appropriate values (because IQTDat::insert_itensor
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
    Real ur = 0.0;
    for(const_iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        { ur += jj->uniqueReal(); }
    return ur;
    }

int IQTensor::
num_index(IndexType t) const
    { 
    int count = 0;
    for(const_iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
        { if(jj->type() == t) ++count; }
    return count;
    }

Real IQTensor::
norm() const
    {
    Real res = 0;
    for(const_iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        { res += sqr(it->norm()); }
    return sqrt(res);
    }

Real IQTensor::
sumels() const
    {
    Real res = 0;
    for(const_iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        { res += it->sumels(); }
    return res;
    }

void IQTensor::
scaleOutNorm() const
    {
    Real f = norm();
    for(const_iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        it->scaleTo(f);
    }

void IQTensor::
scaleTo(LogNumber newscale) const
    {
    for(const_iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        it->scaleTo(newscale);
    }

void IQTensor::
clean(Real min_norm)
    { solo(); p->clean(min_norm); }

void IQTensor::
addindex1(const IQIndex& I)
	{
	if(I.m() != 1) 
	    Error("IQTensor::operator*=(IQIndex): IQIndex must have m == 1.");    
	solo(); 
	p->uninit_rmap();
    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        { it->addindex1(I.index(1)); }
	p->iqindex_.push_back(I);
	}

void IQTensor::
tieIndices(const boost::array<IQIndex,NMAX+1>& indices, int niqind, 
           const IQIndex& tied)
    {
    if(niqind < 1) Error("No IQIndices to tie");

    const int nindex = indices[1].nindex();

    boost::intrusive_ptr<IQTDat> np = new IQTDat();
    np->iqindex_.reserve(p->iqindex_.size()+1-niqind);
    np->iqindex_.push_back(tied);

    int nmatched = 0;
    for(size_t k = 0; k < p->iqindex_.size(); ++k)
        {
        const IQIndex& K = p->iqindex_[k];
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
            np->iqindex_.push_back(K);
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

        for(const_iten_it it=p->itensor.begin(); it != p->itensor.end(); ++it)
            {
            if(!it->hasAllIndex(totie,niqind)) continue;
            np->itensor.push_back(*it);
            np->itensor.back().tieIndices(totie,niqind,tied.index(i));
            }
        }
    np->uninit_rmap();
    p.swap(np);
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

    vector<IQIndex> new_iqinds_;
    new_iqinds_.reserve(p->iqindex_.size()-niqind);

    int nmatched = 0;
    for(size_t k = 0; k < p->iqindex_.size(); ++k)
        {
        const IQIndex& K = p->iqindex_[k];
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
            new_iqinds_.push_back(K);
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

    IQTensor res(new_iqinds_);

    boost::array<Index,NMAX+1> totrace;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 1; n <= niqind; ++n)
            totrace[n] = indices[n].index(i);

        for(const_iten_it it=p->itensor.begin(); it != p->itensor.end(); ++it)
            {
            if(!it->hasAllIndex(totrace,niqind)) continue;
            ITensor tt(*it);
            tt.trace(totrace,niqind);
            res += tt;
            }
        }
    this->operator=(res);
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

int IQTensor::
vecSize() const
	{
    if(this->isNull()) return 0;
	int s = 0;
	for(const_iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    s += jj->vecSize();
	return s;
	}

int IQTensor::
maxSize() const
	{
    int ms = 1;
	for(size_t j = 0; j < p->iqindex_.size(); ++j)
	    ms *= p->iqindex_[j].m();
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
	solo();
	if(vecSize() != v.Length())
	    Error("bad size");
	int off = 1;
	for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    {
	    int d = jj->vecSize();
	    jj->assignFromVec(v.SubVector(off,off+d-1));
	    off += d;
	    }
	}

void IQTensor::
Randomize() 
	{ 
	solo(); 
    for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
        { it->Randomize(); }
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
	for(size_t j = 0; j < p->iqindex_.size(); ++j)
	    cout << p->iqindex_[j] << "\n\n";
	cout << "---------------------------\n" << endl;
	}

void IQTensor::
assignFrom(const IQTensor& other)
	{
    //TODO: account for fermion sign here
    solo();
    if(fabs(uniqueReal()-other.uniqueReal()) > 1E-10)
        {
        PrintIndices((*this));
        PrintIndices(other);
        Error("Can't assign from IQTensor, different IQIndices");
        }

	for(const_iten_it i = other.p->itensor.begin(); i != other.p->itensor.end(); ++i)
	    {
	    ApproxReal r(i->uniqueReal());
	    if(p->has_itensor(r))
            {
            p->rmap[r]->assignFrom(*i);
            }
        else
            {
            p->insert_itensor(r,*i);
            }
	    }
	}

void IQTensor::
conj()
    {
    solo();
    if(!is_complex())
        {
        for(iqind_it jj = p->iqindex_.begin(); jj != p->iqindex_.end(); ++jj)
            { jj->conj(); }
        }
    else
        {
        IQTensor r,i;
        SplitReIm(r,i);
        for(iqind_it jj = r.p->iqindex_.begin(); jj != r.p->iqindex_.end(); ++jj)
            { jj->conj(); }
        for(iqind_it jj = i.p->iqindex_.begin(); jj != i.p->iqindex_.end(); ++jj)
            { jj->conj(); }
        i *= -1.0;
        *this = r * IQTensor::Complex_1() + IQTensor::Complex_i() * i;
        }
    }

void IQTensor::
conj(const IQIndex& I)
    {
    solo();
    Foreach(IQIndex& J, p->iqindex_)
        { 
        if(J == I) J.conj(); 
        }
    }

std::ostream& 
operator<<(std::ostream & s, const IQTensor &t)
    {
    s << "\n----- IQTensor -----\n";
    if(t.isNull())
        {
        s << "(IQTensor is null)\n\n";
        return s;
        }
    s << "IQIndices:\n";
    for(std::vector<IQIndex>::const_iterator
        jj = t.p->iqindex_.begin(); jj != t.p->iqindex_.end(); ++jj)
        { s << "  " << *jj << std::endl; }
    s << "ITensors:\n";
    for(std::list<ITensor>::const_iterator
        it = t.p->itensor.begin(); it != t.p->itensor.end(); ++it)
        { s << "  " << *it << std::endl; }
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
    remove_copy_if(p->iqindex_.begin(),p->iqindex_.end(),std::back_inserter(newreinds),
		    bind2nd(std::equal_to<IQIndex>(),IQIndex::IndReIm()));
    re = IQTensor(newreinds);
    im = re;
    ITensor a,b;
    for(const_iten_it i = p->itensor.begin(); i != p->itensor.end(); ++i)
        {
        i->SplitReIm(a,b);
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
    p->uninit_rmap();

    set<ApproxReal> common_inds;
    
    //Load iqindex_ with those IQIndex's *not* common to *this and other
    static vector<IQIndex> riqind_holder(1000);
    riqind_holder.resize(0);

    for(size_t i = 0; i < p->iqindex_.size(); ++i)
        {
        const IQIndex& I = p->iqindex_[i];
        const_iqind_it f = find(other.p->iqindex_.begin(),other.p->iqindex_.end(),I);
        if(f != other.p->iqindex_.end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() == I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    this->printIndices("*this");
                    other.printIndices("other");
                    cerr << "IQIndex from *this = " << I << endl;
                    cerr << "IQIndex from other = " << *f << endl;
                    Error("Incompatible arrow directions in IQTensor::operator*=.");
                    }
            for(size_t n = 0; n < I.iq().size(); ++n) 
                { common_inds.insert(ApproxReal(I.iq()[n].index.uniqueReal())); }

            common_inds.insert(ApproxReal(I.uniqueReal()));
            }
        else { riqind_holder.push_back(I); }
        }

    for(size_t i = 0; i < other.p->iqindex_.size(); ++i)
    if(!common_inds.count(ApproxReal(other.p->iqindex_[i].uniqueReal())))
        { riqind_holder.push_back(other.p->iqindex_[i]); }

    if(riqind_holder.size() > 1000) cerr << "\nWARNING: in IQTensor::operator*=, riqind_holder had to reallocate.\n\n";
    p->iqindex_.swap(riqind_holder);

    set<ApproxReal> keys;

    list<ITensor> old_itensor; p->itensor.swap(old_itensor);

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
            operator+=(tt);
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

    solo();
    p->uninit_rmap();

    set<ApproxReal> common_inds;
    
    static vector<IQIndex> riqind_holder(1000);
    riqind_holder.resize(0);

    for(size_t i = 0; i < p->iqindex_.size(); ++i)
        {
        const IQIndex& I = p->iqindex_[i];
        const_iqind_it f = find(other.p->iqindex_.begin(),other.p->iqindex_.end(),I);
        if(f != other.p->iqindex_.end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() != I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    this->printIndices("*this");
                    other.printIndices("other");
                    cerr << "IQIndex from *this = " << I << endl;
                    cerr << "IQIndex from other = " << *f << endl;
                    Error("Incompatible arrow directions in IQTensor::operator/=.");
                    }
            for(size_t n = 0; n < I.iq().size(); ++n) 
                { common_inds.insert(ApproxReal(I.iq()[n].index.uniqueReal())); }

            common_inds.insert(ApproxReal(I.uniqueReal()));
            }
        riqind_holder.push_back(I);
        }

    for(size_t i = 0; i < other.p->iqindex_.size(); ++i)
    if(!common_inds.count(ApproxReal(other.p->iqindex_[i].uniqueReal())))
        { riqind_holder.push_back(other.p->iqindex_[i]); }

    if(riqind_holder.size() > 1000) cerr << "\nWARNING: in IQTensor::operator/=, riqind_holder had to reallocate.\n\n";
    p->iqindex_.swap(riqind_holder);

    set<ApproxReal> keys;

    list<ITensor> old_itensor; p->itensor.swap(old_itensor);

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
            operator+=(tt);
            }
        }

    return *this;

    } //IQTensor& IQTensor::operator/=(const IQTensor& other)

//Extracts the real and imaginary parts of the 
//component of a rank 0 tensor (scalar)
void IQTensor::
GetSingComplex(Real& re, Real& im) const
    {
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
    for(const_iqind_it jj = tre.p->iqindex_.begin(); jj != tre.p->iqindex_.end(); ++jj)
        {
        if(*jj != IQTensor::Sing().p->iqindex_[0])
            {
            cout << *this;
            cout << tre;
            Error("bad tre size");
            }
        }
    for(const_iqind_it jj = tim.p->iqindex_.begin(); jj != tim.p->iqindex_.end(); ++jj)
        if(*jj != IQTensor::Sing().p->iqindex_[0])
        { Error("bad tim size"); }

    if(tre.iten_size() == 0)
        { re = 0.0; }
    else
        {
        const ITensor& t = tre.p->itensor.front();
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
        const ITensor& t = tim.p->itensor.front();
        if(t.vecSize() != 1) Error("bad tim dat size");
        im = t.val0();
        }
    }

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    //TODO: account for fermion sign here
    solo(); p->uninit_rmap();
    if(this == &other) 
        {
        for(iten_it it = p->itensor.begin(); it != p->itensor.end(); ++it)
            { *it *= 2; return *this; }
        }

    if(p->iqindex_.size() == 0)		// Automatic initializing a summed IQTensor in a loop
        { return (*this = other); }
    bool complex_this = hasindex(IQIndex::IndReIm()); 
    bool complex_other = other.hasindex(IQIndex::IndReIm()); 
    IQTensor& This(*this);
    if(!complex_this && complex_other)
        return (This = (This * IQTensor::Complex_1()) + other);
    if(complex_this && !complex_other)
        return (This += other * IQTensor::Complex_1());
    if(fabs(This.uniqueReal()-other.uniqueReal()) > 1.0e-11) 
        {
        cout << "This is " << This;
        cout << "other is " << other;
        Error("bad match unique real in IQTensor::operator+=");
        }
    for(const_iten_it it = other.p->itensor.begin(); it != other.p->itensor.end(); ++it)
        { operator+=(*it); }
    return *this;
    }

//
//Automatically convert this IQTensor
//to an ITensor
//
IQTensor::
operator ITensor() const
    {
    //Resulting ITensor's indices are 
    //the Index versions of this's IQIndices
    vector<Index> indices;
    for(size_t j = 0; j < p->iqindex_.size(); ++j)
        {
        indices.push_back(Index(p->iqindex_[j]));
        }
    ITensor res(indices);

    //Loop over ITensors (blocks) within this IQTensor
    for(const_iten_it it = p->itensor.begin();
        it != p->itensor.end();
        ++it)
        {
        ITensor exp(*it);
        //Loop over Index's of the k'th ITensor
        for(int j = 1; j <= it->r(); ++j)
            {
            //Want to transform 'small' into the 
            //Index version of the IQIndex that contains
            //it, with the appropriate offset
            const Index& small = it->index(j);
            //Find the IQIndex that contains 'small'
            const IQIndex* big = 0;
            int offset = -1;
            for(size_t q = 0; q < p->iqindex_.size(); ++q)
                if(p->iqindex_[q].hasindex(small))
                    {
                    big = &(p->iqindex_[q]);
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
solo()
	{
	assert(p != 0);
	if(p->count() != 1)
	    {
	    boost::intrusive_ptr<IQTDat> new_p(new IQTDat(*p));
	    p.swap(new_p);
	    }
	}

Real 
ReSingVal(const IQTensor& x)
    {
    Real re, im;
    x.GetSingComplex(re,im);
    return re;
    }

Real 
Dot(const IQTensor& x, const IQTensor& y, bool doconj)
    {
    IQTensor res(IQTensor::Sing()*(doconj ? conj(x) : x)*y);
    return ReSingVal(res);
    }

void 
Dot(const IQTensor& x, const IQTensor& y, Real& re, Real& im, bool doconj)
    {
    IQTensor res(IQTensor::Sing()*(doconj ? conj(x) : x)*y);
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
