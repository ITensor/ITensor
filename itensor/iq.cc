#include "iq.h"


void DoPrimer::operator()(IQIndex &iqi) const { iqi.doprime(pt,inc); }
void MapPrimer::operator()(IQIndex &iqi) const { iqi.mapprime(plevold,plevnew,pt); }

IQIndexVal IQIndex::operator()(int n) const { return IQIndexVal(*this,n); }

void IQTensor::SplitReIm(IQTensor& re, IQTensor& im) const
{
    if(!hasindex(IQIndReIm))
	{
	IQTensor cop(*this);
	re = cop;
	im = cop;
	im *= 0.0;
	return;
	}
    re = IQTensor();
    remove_copy_if(p->iqindex_.begin(),p->iqindex_.end(),std::back_inserter(re.p->iqindex_),
		    bind2nd(std::equal_to<IQIndex>(),IQIndReIm));
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
    if(hasindex(IQIndReIm) && other.hasindex(IQIndReIm) && !other.hasindex(IQIndReImP)
	    && !other.hasindex(IQIndReImPP))
	{
        static ITensor primer(IndReIm,IndReImP,1.0);
        static ITensor primerP(IndReIm,IndReImPP,1.0);
        static ITensor prod(IndReIm,IndReImP,IndReImPP);
        static IQTensor iqprimer(IQIndReIm,IQIndReImP);
        static IQTensor iqprimerP(IQIndReIm,IQIndReImPP);
        static IQTensor iqprod(IQIndReIm,IQIndReImP,IQIndReImPP);
        static int first = 1;
        if(first)
        {
            IndexVal iv0(IndReIm,1), iv1(IndReImP,1), iv2(IndReImPP,1);
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

    //Handle virtual index
    if(other.viqindex != IQEmptyV)
    if(viqindex == IQEmptyV) viqindex = other.viqindex;
    else
    {
        QN newq = viqindex.dir()*(viqindex.dir()*viqindex.qn(1)+other.viqindex.dir()*other.viqindex.qn(1));
        viqindex = IQIndex((Index) viqindex,viqindex.index(1),newq);
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
            if(f->dir() == I.dir() && f->type() != ReIm && I.type() != ReIm)
            {
                Print(*this);
                Print(other);
                cerr << "IQIndex from *this = " << I << endl;
                cerr << "IQIndex from other = " << *f << endl;
                Error("Incompatible arrow directions in IQTensor::operator*.");
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

    if(riqind_holder.size() > 1000) cerr << "\nWARNING: in IQTensor::operator* riqind_holder had to reallocate.\n\n";
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

} //IQTensor& IQTensor::operator*=(const IQTensor& other)

//Extracts the real and imaginary parts of the 
//component of a rank 0 tensor (scalar)
void IQTensor::GetSingComplex(Real& re, Real& im) const
{
    IQTensor tre,tim;
    SplitReIm(tre,tim);

    //Only IQIndex should be IQIndReIm
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
        if(I.type() != Virtual && I != IQTSing.p->iqindex_[0])
        {
            cout << *this;
            cout << tre;
            Error("bad tre size");
        }
    }
    foreach(const IQIndex& I, tim.p->iqindex_)
    if(I.type() != Virtual && I != IQTSing.p->iqindex_[0])
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
        re = (t.neg() ? -1 : 1) * t.val0() * exp(t.logfac());
	}
    if(tim.iten_size() == 0)
	{ im = 0.0; }
    else
	{
        const ITensor& t = tim.p->itensor.front();
        if(t.vec_size() != 1) Error("bad tim dat size");
        im = (t.neg() ? -1 : 1) * t.val0() * exp(t.logfac());
	}
}

IQTensor& IQTensor::operator+=(const IQTensor& other)
{
    solo(); p->uninit_rmap();
    if(this == &other)
    foreach(ITensor& t, p->itensor) { t *= 2; return *this; }

    if(p->iqindex_.size() == 0)		// Automatic initializing a summed IQTensor in a loop
    { return (*this = other); }
    bool complex_this = hasindex(IQIndReIm); 
    bool complex_other = other.hasindex(IQIndReIm); 
    IQTensor& This(*this);
    if(!complex_this && complex_other)
        return (This = (This * IQComplex_1) + other);
    if(complex_this && !complex_other)
        return (This += other * IQComplex_1);
    if(fabs(This.unique_Real()-other.unique_Real()) > 1.0e-11) 
	{
        cout << "This is " << This;
        cout << "other is " << other;
        Error("bad match unique real in IQTensor::operator+=");
	}
    foreach(const ITensor& t, other.p->itensor) operator+=(t);
    return *this;
}

IQIndex index_in_common(const IQTensor& A, const IQTensor& B, IndexType t)
{
    foreach(const IQIndex& I, A.iqinds())
    if(I.type() == t) if(B.hasindex(I)) return I;

    return IQIndex();
}

