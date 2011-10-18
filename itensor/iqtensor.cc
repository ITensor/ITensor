#include "iqtensor.h"
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::list;
using std::set;
using std::map;
using std::multimap;
using std::ostream;
using std::istream;
using std::pair;
using std::make_pair;

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
	foreach(ITensor& t, itensor) 
	    t.read(s);

	s.read((char*) &size,sizeof(size));
	iqindex_.resize(size);
	foreach(IQIndex& I, iqindex_) 
	    I.read(s);
	}

void IQTDat::
write(ostream& s) const
	{
	size_t size = itensor.size();
	s.write((char*) &size,sizeof(size));
	foreach(const ITensor& t, itensor) 
	    t.write(s);

	size = iqindex_.size();
	s.write((char*) &size,sizeof(size));
	foreach(const IQIndex& I, iqindex_) 
	    I.write(s);
	}

void IQTDat::
init_rmap() const
	{
	if(rmap_init) return;
	for(list<ITensor>::iterator it = itensor.begin(); 
        it != itensor.end(); 
        ++it)
	    rmap[ApproxReal(it->unique_Real())] = it;
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
foreach(const ITensor& t, itensor)
    {
        if(t.norm() >= min_norm)
            nitensor.push_back(t);
    }
itensor.swap(nitensor);
}

void DoPrimer::operator()(IQIndex &iqi) const { iqi.doprime(pt,inc); }
void MapPrimer::operator()(IQIndex &iqi) const { iqi.mapprime(plevold,plevnew,pt); }

void IQTensor::SplitReIm(IQTensor& re, IQTensor& im) const
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

IQTensor& IQTensor::operator*=(const IQTensor& other)
{
    if(this == &other)
        {
        IQTensor cp_oth(other);
        return operator*=(cp_oth);
        }

    if(this->is_null()) 
        Error("'This' IQTensor null in product");

    if(other.is_null()) 
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
        vector<IQIndex>::const_iterator f = find(other.p->iqindex_.begin(),other.p->iqindex_.end(),I);
        if(f != other.p->iqindex_.end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() == I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    this->printIQInds("*this");
                    other.printIQInds("other");
                    cerr << "IQIndex from *this = " << I << endl;
                    cerr << "IQIndex from other = " << *f << endl;
                    Error("Incompatible arrow directions in IQTensor::operator*=.");
                    }
            for(size_t n = 0; n < I.iq().size(); ++n) 
                { common_inds.insert(ApproxReal(I.iq()[n].index.unique_Real())); }

            common_inds.insert(ApproxReal(I.unique_Real()));
            }
        else { riqind_holder.push_back(I); }
        }

    for(size_t i = 0; i < other.p->iqindex_.size(); ++i)
    if(!common_inds.count(ApproxReal(other.p->iqindex_[i].unique_Real())))
        { riqind_holder.push_back(other.p->iqindex_[i]); }

    if(riqind_holder.size() > 1000) cerr << "\nWARNING: in IQTensor::operator*=, riqind_holder had to reallocate.\n\n";
    p->iqindex_.swap(riqind_holder);

    set<ApproxReal> keys;

    list<ITensor> old_itensor; p->itensor.swap(old_itensor);

    //com_this maps the unique_Real of a set of Index's to be contracted over together
    //to those ITensors in *this.itensor having all Index's in that set
    multimap<ApproxReal,const_iten_it> com_this;
    for(const_iten_it tt = old_itensor.begin(); tt != old_itensor.end(); ++tt)
        {
        Real r = 0.0;
        for(int a = 1; a <= tt->r(); ++a)
            {
            if(common_inds.count(ApproxReal(tt->index(a).unique_Real())))
                { r += tt->index(a).unique_Real(); }
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
            if(common_inds.count(ApproxReal(ot->index(b).unique_Real())))
                { r += ot->index(b).unique_Real(); }
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

IQTensor& IQTensor::operator/=(const IQTensor& other)
{
    if(this == &other)
        {
        IQTensor cp_oth(other);
        return operator/=(cp_oth);
        }

    if(this->is_null()) 
        Error("'This' IQTensor null in product");

    if(other.is_null()) 
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
        vector<IQIndex>::const_iterator f = find(other.p->iqindex_.begin(),other.p->iqindex_.end(),I);
        if(f != other.p->iqindex_.end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Globals::checkArrows())
                if(f->dir() != I.dir() && f->type() != ReIm && I.type() != ReIm)
                    {
                    this->printIQInds("*this");
                    other.printIQInds("other");
                    cerr << "IQIndex from *this = " << I << endl;
                    cerr << "IQIndex from other = " << *f << endl;
                    Error("Incompatible arrow directions in IQTensor::operator/=.");
                    }
            for(size_t n = 0; n < I.iq().size(); ++n) 
                { common_inds.insert(ApproxReal(I.iq()[n].index.unique_Real())); }

            common_inds.insert(ApproxReal(I.unique_Real()));
            }
        riqind_holder.push_back(I);
        }

    for(size_t i = 0; i < other.p->iqindex_.size(); ++i)
    if(!common_inds.count(ApproxReal(other.p->iqindex_[i].unique_Real())))
        { riqind_holder.push_back(other.p->iqindex_[i]); }

    if(riqind_holder.size() > 1000) cerr << "\nWARNING: in IQTensor::operator/=, riqind_holder had to reallocate.\n\n";
    p->iqindex_.swap(riqind_holder);

    set<ApproxReal> keys;

    list<ITensor> old_itensor; p->itensor.swap(old_itensor);

    //com_this maps the unique_Real of a set of Index's to be summed over together
    //to those ITensors in *this.itensor having all Index's in that set
    multimap<ApproxReal,const_iten_it> com_this;
    for(const_iten_it tt = old_itensor.begin(); tt != old_itensor.end(); ++tt)
        {
        Real r = 0.0;
        for(int a = 1; a <= tt->r(); ++a)
            {
            Real ur = tt->index(a).unique_Real();
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
            Real ur = ot->index(b).unique_Real();
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
void IQTensor::GetSingComplex(Real& re, Real& im) const
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
    foreach(const IQIndex& I, tre.p->iqindex_)
    {
        if(I != IQTensor::Sing().p->iqindex_[0])
        {
            cout << *this;
            cout << tre;
            Error("bad tre size");
        }
    }
    foreach(const IQIndex& I, tim.p->iqindex_)
    if(I != IQTensor::Sing().p->iqindex_[0])
    { Error("bad tim size"); }

    if(tre.iten_size() == 0)
    { re = 0.0; }
    else
	{
        const ITensor& t = tre.p->itensor.front();
        if(t.vec_size() != 1) 
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
        if(t.vec_size() != 1) Error("bad tim dat size");
        im = t.val0();
	}
}

IQTensor& IQTensor::operator+=(const IQTensor& other)
{
    solo(); p->uninit_rmap();
    if(this == &other) {
        foreach(ITensor& t, p->itensor) { t *= 2; return *this; }
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
    if(fabs(This.unique_Real()-other.unique_Real()) > 1.0e-11) 
	{
        cout << "This is " << This;
        cout << "other is " << other;
        Error("bad match unique real in IQTensor::operator+=");
	}
    foreach(const ITensor& t, other.p->itensor) operator+=(t);
    return *this;
}

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
    for(list<ITensor>::const_iterator it = p->itensor.begin();
        it != p->itensor.end();
        ++it)
        {
        ITensor exp(*it),nexp;
        //Loop over Index's of the k'th ITensor
        for(int j = 1; j <= it->r(); ++j)
            {
            //Want to transform 'small' into the 
            //Index version of the IQIndex that contains
            //it, with the appropriate offset
            const Index& small = it->index(j);
            //Find the IQIndex that contains 'small'
            const IQIndex* big;
            int offset = -1;
            for(size_t q = 0; q < p->iqindex_.size(); ++q)
                if(p->iqindex_[q].hasindex(small))
                    {
                    big = &(p->iqindex_[q]);
                    offset = big->offset(small);
                    break;
                    }
            exp.expandIndex(small,*big,offset,nexp);
            exp = nexp;
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
