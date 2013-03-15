//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqtsparse.h"
#include <set>
using namespace std;
using boost::format;
using boost::array;
using boost::shared_ptr;
using boost::make_shared;

//
// IQTSDat
//

IQTSDat::
IQTSDat()
    :
    init(false)
    {
    }

void IQTSDat::
makeCopyOf(const IQTSDat& other)
    {
    init = false;
    its_ = other.its_;
    }

IQTSDat::
IQTSDat(istream& s)
    :
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
    rmap.clear();
    init = false;
    }

void IQTSDat::
scaleTo(const LogNumber& newscale) const
    {
    Foreach(ITSparse& s, its_)
        s.scaleTo(newscale);
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

const shared_ptr<IQTSDat>& IQTSDat::
Null()
    {
    static shared_ptr<IQTSDat> Null_ = make_shared<IQTSDat>();
    return Null_;
    }

//
// IQTSparse
//

IQTSparse::
IQTSparse()
    :
    is_(IndexSet<IQIndex>::Null()),
    d_(IQTSDat::Null())
    { }

IQTSparse::
IQTSparse(const IQIndex& i1)
    :
    is_(make_shared<IndexSet<IQIndex> >(i1)),
    d_(make_shared<IQTSDat>())
    { 
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2)
    :
    is_(make_shared<IndexSet<IQIndex> >(i1,i2)),
    d_(make_shared<IQTSDat>())
    { 
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2,
          Real r)
    :
    is_(make_shared<IndexSet<IQIndex> >(i1,i2)),
    d_(make_shared<IQTSDat>())
    { 
    if(i1.nindex() != i2.nindex())
        {
        Print(i1);
        Print(i2);
        Error("IQIndex's have different number of qn blocks");
        }
    for(int j = 1; j <= i1.nindex(); ++j)
        operator+=(ITSparse(i1.index(j),i2.index(j),r));
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3)
    :
    is_(make_shared<IndexSet<IQIndex> >(i1,i2,i3)),
    d_(make_shared<IQTSDat>())
    { 
    }

bool IQTSparse::
isNull() const { return (d_ == IQTSDat::Null() || is_ == IndexSet<IQIndex>::Null()); }

bool IQTSparse::
isNotNull() const { return !(d_ == IQTSDat::Null() || is_ == IndexSet<IQIndex>::Null()); }


IQTSparse& IQTSparse::
operator+=(const ITSparse& s)
    {
    ncblocks().insert_add(s);
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

    Foreach(const ITSparse& s, other.blocks())
        {
        ncblocks().insert_add(s);
        }

    return *this;
    }

IQTSparse& IQTSparse::
operator*=(Real fac)
    { 
    soloDat();

    if(fac == 0)
        {
        ncblocks().clear();
        return *this;
        }

    Foreach(ITSparse& s, ncblocks())
        {
        s *= fac;
        }
    return *this;
    }

IQTSparse& IQTSparse::
operator*=(const LogNumber& fac)
    { 
    soloDat();

    if(fac == 0)
        {
        ncblocks().clear();
        return *this;
        }

    Foreach(ITSparse& s, ncblocks())
        {
        s *= fac;
        }
    return *this;
    }

/*
Real& IQTSparse::
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

    if(!blocks().has_itensor(r))
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

    return (ncblocks().get(r)).operator()(iv1.blockIndexVal(),
                                       iv2.blockIndexVal(),
                                       iv3.blockIndexVal(),
                                       iv4.blockIndexVal(),
                                       iv5.blockIndexVal(),
                                       iv6.blockIndexVal(),
                                       iv7.blockIndexVal(),
                                       iv8.blockIndexVal());
	}
    */

//
// Primelevel Methods 
//

void IQTSparse::
noprime(IndexType type)
    { 
    solo();

    is_->noprime(type); 

    Foreach(ITSparse& t, ncblocks())
        {
        t.noprime(type);
        }
    }

void IQTSparse::
prime(IndexType type, int inc) 
    { 
    solo();

    is_->prime(type,inc);

    Foreach(ITSparse& t, ncblocks())
        {
        t.prime(type,inc);
        }
    }

void IQTSparse::
mapprime(int plevold, int plevnew, IndexType type)
    { 
    solo();

    is_->mapprime(plevold,plevnew,type); 

    Foreach(ITSparse& t, ncblocks())
        {
        t.mapprime(plevold,plevnew,type);
        }
    }

void IQTSparse::
prime(const IQIndex& I, int inc)
    {
    solo();

    is_->prime(I,inc);

    Foreach(ITSparse& t, ncblocks())
    Foreach(const Index& i, I.indices())
        {
        if(hasindex(t,i))
            t.prime(i,inc);
        }
    }

void IQTSparse::
noprime(const IQIndex& I) 
    { 
    solo();

    is_->noprime(I); 

    Foreach(ITSparse& t, ncblocks())
    Foreach(const Index& i, I.indices())
        {
        if(hasindex(t,i))
            t.noprime(i);
        }
    }

void IQTSparse::
conj()
    {
    if(!isComplex(*this))
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

void IQTSparse::
pseudoInvert(Real cutoff)
    {
    soloDat();

    Foreach(ITSparse& s, ncblocks())
        { 
        s.pseudoInvert(cutoff);
        }
    }

Real IQTSparse::
norm() const
    {
    Real res = 0;
    Foreach(const ITSparse& s, blocks())
        {
        res += sqr(s.norm());
        }
    return sqrt(res);
    }

void IQTSparse::
scaleOutNorm() const
	{
    Real nrm = norm();
    blocks().scaleTo(nrm);
	}

void IQTSparse::
scaleTo(const LogNumber& newscale) const
	{
    blocks().scaleTo(newscale);
	}

void IQTSparse::
print(std::string name,Printdat pdat) const 
	{ 
    bool savep = Global::printdat();
    Global::printdat() = (pdat==ShowData); 
	std::cerr << "\n" << name << " =\n" << *this << "\n"; 
    Global::printdat() = savep;
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
    is_ = make_shared<IndexSet<IQIndex> >();
    is_->read(s);
    d_ = make_shared<IQTSDat>();
    d_->read(s);
    }

void IQTSparse::
write(std::ostream& s) const
    {
	bool null_ = isNull();
	s.write((char*) &null_,sizeof(null_));
	if(null_) return;
    is_->write(s);
	blocks().write(s);
    }

void IQTSparse::
soloDat()
    {
	if(d_ == 0)
        {
        Error("IQTSparse is null");
        }

    if(!d_.unique())
	    {
        const IQTSDat& old_dat(*d_);
        d_ = make_shared<IQTSDat>();
        d_->makeCopyOf(old_dat);
	    }
    }

void IQTSparse::
soloIndex()
    {
	if(!is_)
        Error("IQTSparse is null");

	if(!is_.unique())
        {
        const IndexSet<IQIndex>& old_is(*is_);
        is_ = make_shared<IndexSet<IQIndex> >();
        *is_ = old_is;
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

    if(hasindex(S,IQIndex::IndReIm()) && hasindex(T,IQIndex::IndReIm()) && !hasindex(T,IQIndex::IndReImP())
	    && !hasindex(T,IQIndex::IndReImPP()) && !hasindex(S,IQIndex::IndReImP()) && !hasindex(S,IQIndex::IndReImPP()))
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
        IQTensor::const_iqind_it f = find(T.is_->begin(),T.is_->end(),I);
        if(f != T.is_->end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Global::checkArrows())
                if(f->dir() == I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    Print(S.indices());
                    Print(T.indices());
                    cerr << "IQIndex from S = " << I << endl;
                    cerr << "IQIndex from T = " << *f << endl;
                    cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                    throw ArrowError("Incompatible arrow directions in IQTensor::operator*=.");
                    }
            Foreach(const Index& i, I.indices())
                { common_inds.insert(ApproxReal(i.uniqueReal())); }

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

    IQTDat::StorageT old_itensor; 
    res.dat.nc().swap(old_itensor);

    multimap<ApproxReal,IQTSDat::const_iterator> com_S;
    for(IQTSDat::const_iterator tt = S.blocks().begin(); tt != S.blocks().end(); ++tt)
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

    multimap<ApproxReal,IQTDat::const_iterator> com_T;
    for(IQTDat::const_iterator ot = T.blocks().begin(); ot != T.blocks().end(); ++ot)
        {
        Real r = 0.0;
        Foreach(const Index& I, ot->indices())
            {
            if(common_inds.count(ApproxReal(I.uniqueReal())))
                { r += I.uniqueReal(); }
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
            if(tt.scale().sign() != 0)
                res.dat.nc().insert_add(tt);
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
    Foreach(const ITSparse& t, T.blocks())
        { s << "  " << t << std::endl; }
    s << "-------------------" << "\n\n";
    return s;
    }
