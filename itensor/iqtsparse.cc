//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqtsparse.h"
#include <set>
using namespace std;
using boost::format;
using boost::array;

//
// IQTSDat
//

IQTSDat::
IQTSDat()
    :
    numref(0),
    init(false)
    {
    }

IQTSDat::
IQTSDat(istream& s)
    :
    numref(0),
    init(false)
    {
    read(s);
    }

void IQTSDat::
insert_add(const ITSparse& s)
    {
    init_rmap();

    ApproxReal r(s.uniqueReal());
    if(rmap.count(r) == 1)
        {
        *rmap[r] += s;
        return;
        }
    else
        {
        its_.push_front(s);
        rmap[r] = its_.begin();
        }
    }

void IQTSDat::
clear()
    {
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("clear called on shared IQTSDat");
        }
#endif
    rmap.clear();
    its_.clear();
    }

void IQTSDat::
init_rmap() const
    {
    if(init) return;

    for(iterator it = its_.begin(); it != its_.end(); ++it)
        rmap[ApproxReal(it->uniqueReal())] = it;

    init = true;
    }

void IQTSDat::
uninit_rmap() const
    {
#ifdef DEBUG
    if(numref > 1)
        {
        Print(numref);
        Error("uninit_rmap called on shared IQTSDat");
        }
#endif
    rmap.clear();
    init = false;
    }

void IQTSDat::
read(istream& s)
    {
    uninit_rmap();
	size_t size;
	s.read((char*) &size,sizeof(size));
	its_.resize(size);
    for(iterator it = its_.begin(); it != its_.end(); ++it)
        { it->read(s); }
    }

void IQTSDat::
write(ostream& s) const
    {
	size_t size = its_.size();
	s.write((char*) &size,sizeof(size));
    for(const_iterator it = its_.begin(); it != its_.end(); ++it)
        { it->write(s); }
    }


//
// IQTSparse
//

IQTSparse::
IQTSparse()
    :
    is_(0),
    d_(0)
    { }

IQTSparse::
IQTSparse(const IQIndex& i1)
    :
    is_(new IQIndexSet(i1)),
    d_(new IQTSDat())
    { 
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2)
    :
    is_(new IQIndexSet(i1,i2)),
    d_(new IQTSDat())
    { 
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3)
    :
    is_(new IQIndexSet(i1,i2,i3)),
    d_(new IQTSDat())
    { 
    }


IQTSparse& IQTSparse::
operator+=(const ITSparse& s)
    {
    d_->insert_add(s);
    return *this;
    }

IQTSparse& IQTSparse::
operator+=(const IQTSparse& other)
    {
    if(this == &other)
        {
        operator*=(2);
        return *this;
        }

    Foreach(const ITSparse& s, *(other.d_))
        {
        d_->insert_add(s);
        }

    return *this;
    }

IQTSparse& IQTSparse::
operator*=(Real fac)
    { 
    soloDat();

    if(fac == 0)
        {
        d_->clear();
        return *this;
        }

    Foreach(ITSparse& s, *d_)
        {
        s *= fac;
        }
    return *this;
    }

IQIndex IQTSparse::
findtype(IndexType t) const 
    { 
    if(is_ == 0)
        {
        Error("IQTSparse is null");
        }
    return is_->findtype(t); 
    }

bool IQTSparse::
findtype(IndexType t, IQIndex& I) const 
    { 
    if(is_ == 0) return false;

    return is_->findtype(t,I); 
    }

int IQTSparse::
findindex(const IQIndex& I) const 
    { 
    if(is_ == 0)
        {
        Error("IQTSparse is null");
        }
    return is_->findindex(I); 
    }

bool IQTSparse::
has_common_index(const IQTSparse& other) const
    { 
    if(is_ == 0) return false;

    return is_->has_common_index(*other.is_); 
    }

bool IQTSparse::
hasindex(const IQIndex& I) const 
    { 
    if(is_ == 0) return false;

    return is_->hasindex(I); 
    }

//
// Primelevel Methods 
//

void IQTSparse::
noprime(PrimeType pt)
    { 
    solo();

    is_->noprime(pt); 

    Foreach(ITSparse& t, *(d_))
        {
        t.noprime(pt);
        }
    }

void IQTSparse::
doprime(PrimeType pt, int inc) 
    { 
    solo();

    is_->doprime(pt,inc);

    Foreach(ITSparse& t, *(d_))
        {
        t.doprime(pt,inc);
        }
    }

void IQTSparse::
mapprime(int plevold, int plevnew, PrimeType pt)
    { 
    solo();

    is_->mapprime(plevold,plevnew,pt); 

    Foreach(ITSparse& t, *(d_))
        {
        t.mapprime(plevold,plevnew,pt);
        }
    }

void IQTSparse::
mapprimeind(const IQIndex& I, int plevold, int plevnew, 
            PrimeType pt)
    { 
    solo();

    is_->mapprimeind(I,plevold,plevnew,pt); 

    Foreach(ITSparse& t, *(d_))
        {
        t.mapprimeind(I,plevold,plevnew,pt);
        }
    }

void IQTSparse::
primeind(const IQIndex& I, int inc)
    {
    solo();

    is_->primeind(I,inc);

    Foreach(ITSparse& t, *(d_))
    for(std::vector<inqn>::const_iterator x = I.iq().begin(); 
            x != I.iq().end(); ++x)
        {
        if(t.hasindex(x->index))
            t.primeind(x->index,inc);
        }
    }

void IQTSparse::
noprimeind(const IQIndex& I) 
    { 
    solo();

    is_->noprimeind(I); 

    Foreach(ITSparse& t, *(d_))
    for(std::vector<inqn>::const_iterator x = I.iq().begin(); 
            x != I.iq().end(); ++x)
        {
        if(t.hasindex(x->index))
            t.noprimeind(x->index);
        }
    }

void IQTSparse::
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
        Error("Complex IQTSparse not yet implemented");
        }
    }

Real IQTSparse::
norm() const
    {
    Real res = 0;
    Foreach(const ITSparse& s, (*d_))
        {
        res += sqr(s.norm());
        }
    return sqrt(res);
    }

void IQTSparse::
scaleOutNorm() const
	{
    Real nrm = norm();
    Foreach(ITSparse& s, (*d_))
        {
        s.scaleTo(nrm);
        }
	}

void IQTSparse::
scaleTo(const LogNumber& newscale) const
	{
    Foreach(ITSparse& s, (*d_))
        {
        s.scaleTo(newscale);
        }
	}

void IQTSparse::
print(std::string name,Printdat pdat) const 
	{ 
    bool savep = Globals::printdat();
    Globals::printdat() = (pdat==ShowData); 
	std::cerr << "\n" << name << " =\n" << *this << "\n"; 
    Globals::printdat() = savep;
	}

void IQTSparse::
printIndices(const std::string& name) const
	{ 
	cout << "\n" << name << " (IQIndices only) = \n";
    if(this->isNull())
        {
        cout << "    [IQTSparse is null]" << endl;
        return;
        }
	for(int j = 1; j <= is_->r(); ++j)
	    cout << is_->index(j) << "\n\n";
	cout << "---------------------------\n" << endl;
	}

void IQTSparse::
read(std::istream& s)
    {
    bool null_;
    s.read((char*) &null_,sizeof(null_));
    if(null_) 
        { *this = IQTSparse(); return; }
    is_ = new IQIndexSet(s);
    d_ = new IQTSDat(s);
    }

void IQTSparse::
write(std::ostream& s) const
    {
	bool null_ = isNull();
	s.write((char*) &null_,sizeof(null_));
	if(null_) return;
    is_->write(s);
	d_->write(s);
    }

void IQTSparse::
soloDat()
    {
	if(d_ == 0)
        {
        Error("IQTSparse is null");
        }

	if(d_->count() != 1)
	    {
	    boost::intrusive_ptr<IQTSDat> new_d_(new IQTSDat(*d_));
	    d_.swap(new_d_);
	    }
    }

void IQTSparse::
soloIndex()
    {
	if(is_ == 0)
        {
        Error("IQTSparse is null");
        }

	if(is_->count() != 1)
	    {
	    boost::intrusive_ptr<IQIndexSet> 
            new_is_(new IQIndexSet(*is_));
	    is_.swap(new_is_);
	    }
    }

void IQTSparse::
solo()
    {
    soloIndex();
    soloDat();
    }

void
product(const IQTSparse& S, const IQTensor& T, IQTensor& res)
    {
    if(S.isNull()) 
        Error("IQTSparse null in product");

    if(T.isNull()) 
        Error("Multiplying by null IQTensor");

    if(S.hasindex(IQIndex::IndReIm()) && T.hasindex(IQIndex::IndReIm()) && !T.hasindex(IQIndex::IndReImP())
	    && !T.hasindex(IQIndex::IndReImPP()) && !S.hasindex(IQIndex::IndReImP()) && !S.hasindex(IQIndex::IndReImPP()))
        {
        Error("Complex IQTSparse not yet implemented");
        }

    set<ApproxReal> common_inds;
    
    //Load iqindex_ with those IQIndex's *not* common to *this and other
    static vector<IQIndex> riqind_holder;
    riqind_holder.resize(0);

    for(int i = 1; i <= S.is_->r(); ++i)
        {
        const IQIndex& I = S.is_->index(i);
        IQTensor::const_iqind_it f = find(T.is_->index_.begin(),T.is_->index_.end(),I);
        if(f != T.is_->index_.end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() == I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    PrintIndices(S);
                    PrintIndices(T);
                    cerr << "IQIndex from S = " << I << endl;
                    cerr << "IQIndex from T = " << *f << endl;
                    Error("Incompatible arrow directions in IQTensor::operator*=.");
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

    for(int i = 1; i <= T.is_->r(); ++i)
        {
        const IQIndex& I = T.is_->index(i);
        if(!common_inds.count(ApproxReal(I.uniqueReal())))
            { 
            riqind_holder.push_back(I); 
            }
        }

    res = IQTensor(riqind_holder);

    set<ApproxReal> keys;

    list<ITensor> old_itensor; 
    res.p->uninit_rmap();
    res.p->itensor.swap(old_itensor);

    multimap<ApproxReal,IQTSDat::const_iterator> com_S;
    for(IQTSDat::const_iterator tt = S.d_->begin(); tt != S.d_->end(); ++tt)
        {
        Real r = 0.0;
        for(int a = 1; a <= tt->r(); ++a)
            {
            if(common_inds.count(ApproxReal(tt->index(a).uniqueReal())))
                { r += tt->index(a).uniqueReal(); }
            }
        com_S.insert(make_pair(ApproxReal(r),tt));
        keys.insert(ApproxReal(r));
        }

    multimap<ApproxReal,IQTensor::const_iten_it> com_T;
    for(IQTensor::const_iten_it ot = T.const_iten_begin(); ot != T.const_iten_end(); ++ot)
        {
        Real r = 0.0;
        for(int b = 1; b <= ot->r(); ++b)
            {
            if(common_inds.count(ApproxReal(ot->index(b).uniqueReal())))
                { r += ot->index(b).uniqueReal(); }
            }
        com_T.insert(make_pair(ApproxReal(r),ot));
        keys.insert(ApproxReal(r));
        }

    typedef multimap<ApproxReal,IQTensor::const_iten_it>::iterator 
    rit;
    pair<rit,rit> rrange;

    typedef multimap<ApproxReal,IQTSDat::const_iterator>::iterator 
    lit;
    pair<lit,lit> lrange;

    ITensor tt;
    for(set<ApproxReal>::iterator k = keys.begin(); k != keys.end(); ++k)
        {
        //Equal range returns the begin and end iterators for the sequence
        //corresponding to multimap[key] as a pair
        lrange = com_S.equal_range(*k);
        rrange = com_T.equal_range(*k);

        //Iterate over all ITensors in *this and other sharing
        //the set of contracted Index's corresponding to k
        for(lit ll = lrange.first; ll != lrange.second; ++ll)
        for(rit rr = rrange.first; rr != rrange.second; ++rr)
            {
            //Multiply the ITensors and add into res
            tt = *(ll->second) * *(rr->second);
            if(!tt.scale().isRealZero())
                res.p->insert_add(tt);
            }
        }

    } // void product(const IQTSparse& S, const IQTensor& T, IQTensor& res)

ostream& 
operator<<(ostream & s, const IQTSparse& T)
    {
    s << "\n----- IQTSparse -----\n";
    if(T.isNull())
        {
        s << "(IQTSparse is null)\n\n";
        return s;
        }
    s << "IQIndices:\n";
    for(int k = 1; k <= T.r(); ++k)
        { s << "  " << T.index(k) << std::endl; }
    s << "ITSparse blocks:\n";
    Foreach(const ITSparse& t, (*T.d_))
        { s << "  " << t << std::endl; }
    s << "-------------------" << "\n\n";
    return s;
    }
