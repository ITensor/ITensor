//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"
#include "qcounter.h"
using namespace std;
using boost::format;
using boost::array;
//using boost::shared_ptr;
//using boost::make_shared;

IQTensor::Data::
Data()
    : p(boost::make_shared<IQTDat<ITensor> >())
    { }

IQTensor::Data::
Data(const IQTDatPtr& p_)
    : p(p_)
    { }

void IQTensor::Data::
solo()
    {
#ifdef DEBUG
	if(!p) Error("IQTensor is null");
#endif
	if(!p.unique())
        {
        p = boost::make_shared<IQTDat<ITensor> >(*p);
        }
	}



//
// IQTensor
//

bool IQTensor::
isNull() const 
    { 
    return &dat() == IQTDat<ITensor>::Null().get(); 
    }

bool IQTensor::
isComplex() const
    {
    Foreach(const ITensor& t, dat())
        {
        if(t.isComplex()) return true;
        }
    return false;
    }

int IQTensor::
r() const 
    { 
    return is_->r(); 
    }

/*
const IQIndex& IQTensor::
index(int j) const 
    { 
    if(!is_) Error("IQTensor is null");
    return is_->index(j); 
    }
    */

bool IQTensor::
empty() const { return dat().empty(); }

//----------------------------------------------------
//IQTensor: Constructors 

IQTensor::
IQTensor() 
    : 
    is_(IndexSet<IQIndex>::Null()),
    dat(IQTDat<ITensor>::Null()) 
    { }

IQTensor::
IQTensor(Real val) 
    : 
    is_(boost::make_shared<IndexSet<IQIndex> >())
    { 
    operator+=(ITensor(val));
    }

IQTensor::
IQTensor(const IQIndex& i1) 
    : 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2) 
    : 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) 
	: 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2,i3))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4) 
    : 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2,i3,i4))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5) 
    : 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2,i3,i4,i5))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6)
	: 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2,i3,i4,i5,i6))
	{ 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7)
	: 
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2,i3,i4,i5,i6,i7))
	{ 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7,const IQIndex& i8)
	: 
    is_(new IndexSet<IQIndex>(i1,i2,i3,i4,i5,i6,i7,i8))
	{ 
    }

IQTensor::
IQTensor(std::vector<IQIndex>& iqinds_) 
	: 
    is_(new IndexSet<IQIndex>(iqinds_))
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
    is_(boost::make_shared<IndexSet<IQIndex> >(IQIndex(iv1)))
	{ 
	operator()(iv1) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2) 
	: 
    is_(boost::make_shared<IndexSet<IQIndex> >(iv1,iv2))
	{ 
    operator()(iv1,iv2) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
         const IQIndexVal& iv3)
	: 
    is_(boost::make_shared<IndexSet<IQIndex> >(iv1,iv2,iv3))
	{ 
    operator()(iv1,iv2,iv3) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
         const IQIndexVal& iv3, const IQIndexVal& iv4)
	: 
    is_(boost::make_shared<IndexSet<IQIndex> >(iv1,iv2,iv3,iv4))
	{ 
    operator()(iv1,iv2,iv3,iv4) = 1;
	}

/* Deprecated (Jul 15, 2013 EMS)
IQTensor::
IQTensor(IndexType type,const IQTensor& other) 
    : 
    is_(other.is_),
    dat(other.dat)
    { prime(type); }
*/

void IQTensor::
read(std::istream& s)
    {
    bool null_;
    s.read((char*) &null_,sizeof(null_));
    if(null_) 
        { *this = IQTensor(); return; }
    is_ = boost::make_shared<IndexSet<IQIndex> >();
    is_->read(s);
    dat = Data();
    dat.nc().read(s);
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
    dat.solo();

    if(fac == 0) 
        { 
        dat.nc().clear(); 
        return *this; 
        }

    Foreach(ITensor& t, dat.nc())
        {
        t *= fac;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator/=(Real fac) 
    { 
    dat.solo();

    Foreach(ITensor& t, dat.nc())
        {
        t /= fac;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator*=(Complex z) 
    { 
    dat.solo();

    if(z.real() == 0 && z.imag() == 0) 
        { 
        dat.nc().clear(); 
        return *this; 
        }

    Foreach(ITensor& t, dat.nc())
        {
        t *= z;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator*=(const LogNumber& lgnum) 
    { 
    dat.solo();

    Foreach(ITensor& t, dat.nc())
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
        dat.solo();
        dat.nc().insert(t);
        }
    }

IQTensor& IQTensor::
operator+=(const ITensor& t) 
    { 
#ifdef DEBUG
    if(!this->empty())
        {
        const
        QN d = div(*this);
        QN q;
        Foreach(const Index& i, t.indices())
            {
            q += qn(*this,i)*dir(*this,i);
            }
        if(q != d)
            {
            //Print(d);
            //Print(q);
            throw ITError("New ITensor block has different divergence from IQTensor.");
            }
        }
#endif
    if(t.scale().sign() != 0)
        {
        dat.solo();
        dat.nc().insert_add(t);
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
    dat.solo();
    boost::array<IQIndexVal,NMAX+1> iv 
        = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};

    int nn = 0; 
    IndexSet<Index> is;
    while(GET(iv,nn+1) != IQIndexVal::Null()) 
        {
        ++nn;
        if(!hasindex(*this,iv.at(nn))) 
            Error("IQTensor::operator(): IQIndex not found.");
        is.addindex(iv[nn].indexqn());
        }
    if(nn != r()) 
        {
        Error("Wrong number of IQIndexVals provided");
        }

    return (dat.nc().get(is)).operator()(iv1.blockIndexVal(),
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

    int nn = 0; 
    IndexSet<Index> is;
    while(GET(iv,nn+1) != IQIndexVal::Null()) 
        {
        //ur += GET(iv,++nn).indexqn().uniqueReal(); 
        if(!hasindex(*this,iv.at(nn))) 
            Error("IQTensor::operator(): IQIndex not found.");
        is.addindex(iv[nn].indexqn());
        }
    if(nn != r()) 
        Error("Wrong number of IQIndexVals provided");

    if(!dat().hasBlock(is))
        {
        return 0.;
        }
    else
        {
        return (dat().get(is)).operator()(iv1.blockIndexVal(),
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





IQTensor& IQTensor::
noprime(IndexType type)
	{
	solo();

    is_->noprime(type);

    Foreach(ITensor& t, dat.nc())
        t.noprime(type); 
    return *this;
	} 

IQTensor& IQTensor::
prime(IndexType type, int inc)
	{
	solo();

    is_->prime(type,inc);

    Foreach(ITensor& t, dat.nc())
	    t.prime(type,inc);
    return *this;
	}

IQTensor& IQTensor::
mapprime(int plevold, int plevnew, IndexType type)
    {
    solo();

    is_->mapprime(plevold,plevnew,type);

    Foreach(ITensor& t, dat.nc())
	    t.mapprime(plevold,plevnew,type);
    return *this;
	}

IQTensor& IQTensor::
prime(const IQIndex& I, int inc)
	{
	solo();

    is_->prime(I,inc);

    Foreach(ITensor& t, dat.nc())
    Foreach(const Index& i, I.indices())
        {
		if(hasindex(t,i)) 
		    t.prime(i,inc);
        }
    return *this;
	}

IQTensor& IQTensor::
noprime(const IQIndex& I)
	{
	solo();

    is_->noprime(I);

    Foreach(ITensor& t, dat.nc())
    Foreach(const Index& i, I.indices())
        {
        if(hasindex(t,i)) 
            t.noprime(i);
        }
    return *this;
	}




LogNumber IQTensor::
normLogNum() const
    {
    if(dat().empty()) return LogNumber(0);

    Real maxLogNum = -maxlogdouble;
    Foreach(const ITensor& t, dat())
        { 
        if(t.scale().logNum() > maxLogNum)
            maxLogNum = t.scale().logNum();
        }

    //Add sum of squares of the block norms with exp(2*maxLogNum) scaled out, 
    //using lognorm as a temporary
    Real lognorm = 0;
    Foreach(const ITensor& t, dat())
        { 
        if(t.scale().sign() != 0)
            {
            lognorm += exp(2*(t.scale().logNum()-maxLogNum))*sqr(t.normNoScale()); 
            }
        }
    //Take the sqrt (the 0.5 factor) and log of lognorm, then add maxLogNum back in
    lognorm = maxLogNum + 0.5*log(lognorm);

    return LogNumber(lognorm,+1);
    }

Real IQTensor::
norm() const
    {
    return normLogNum().real0();
    }


Real 
sumels(const IQTensor& T)
    {
    Real sum = 0;
    Foreach(const ITensor& t, T.blocks())
        { 
        sum += sumels(t); 
        }
    return sum;
    }

void IQTensor::
scaleOutNorm()
    {
    dat.solo(); 
    const LogNumber newscale = normLogNum();
    Foreach(ITensor& t, dat.nc())
        t.scaleTo(newscale);
    }

void IQTensor::
scaleTo(const LogNumber& newscale)
    {
    dat.solo(); 
    Foreach(ITensor& t, dat.nc())
        t.scaleTo(newscale);
    }

void IQTensor::
clean(Real min_norm)
    { 
    dat.solo(); 
    dat.nc().clean(min_norm); 
    }


void IQTensor::
tieIndices(const boost::array<IQIndex,NMAX>& indices, int niqind, 
           const IQIndex& tied)
    {
    if(niqind < 1) Error("No IQIndices to tie");

    const int nindex = indices[0].nindex();

    Data ndat;
    boost::shared_ptr<IndexSet<IQIndex> > nis_ = boost::make_shared<IndexSet<IQIndex> >(tied);

    int nmatched = 0;
    for(int k = 1; k <= is_->r(); ++k)
        {
        const IQIndex& K = is_->index(k);
        bool K_is_tied = false;
        for(int j = 0; j < niqind; ++j)
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
        Print(this->indices());
        cout << "Indices to tie = " << endl;
        for(int j = 0; j < niqind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find IQIndex to tie");
        }

    boost::array<Index,NMAX> totie;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 0; n < niqind; ++n)
            totie[n] = indices[n].index(i);

        Foreach(const ITensor& t, dat())
            {
            bool has_all = true;
            for(int n = 0; n < niqind; ++n)
                {
                if(!hasindex(t,totie[n]))
                    {
                    has_all = false;
                    break;
                    }
                }
            if(has_all)
                {
                ITensor nt(t);
                nt.tieIndices(totie,niqind,tied.index(i));
                ndat.nc().insert_add(nt);
                }
            }
        }
    dat.swap(ndat);
    is_.swap(nis_);
    }

void IQTensor::
tieIndices(const IQIndex& i1, const IQIndex& i2, const IQIndex& tied)
    {
    boost::array<IQIndex,NMAX> inds =
        {{ i1, i2, 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null() }};

    tieIndices(inds,2,tied);
    }

IQTensor& IQTensor::
trace(const boost::array<IQIndex,NMAX>& indices, int niqind)
    {
    if(niqind < 0)
        {
        niqind = 0;
        while(indices[niqind] != IQIndex::Null()) ++niqind;
        }

    if(niqind < 1) Error("No IQIndices to trace");

    const int nindex = indices[0].nindex();
    const int tm = indices[0].m();

    Data ndat;
    boost::shared_ptr<IndexSet<IQIndex> > nis_ = boost::make_shared<IndexSet<IQIndex> >();

    int nmatched = 0;
    for(int k = 1; k <= is_->r(); ++k)
        {
        const IQIndex& K = is_->index(k);
        bool K_traced = false;
        for(int j = 0; j < niqind; ++j)
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
        Print(this->indices());
        cout << "Indices to trace = " << endl;
        for(int j = 0; j < niqind; ++j)
            cout << indices[j] << endl;
        Error("Couldn't find IQIndex to trace");
        }

    boost::array<Index,NMAX> totrace;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 0; n < niqind; ++n)
            totrace[n] = indices[n].index(i);

        Foreach(const ITensor& t, dat())
            {
            bool has_all = true;
            for(int n = 0; n < niqind; ++n)
                {
                if(!hasindex(t,totrace[n]))
                    {
                    has_all = false;
                    break;
                    }
                }
            if(has_all)
                {
                ITensor tt(t);
                tt.trace(totrace,niqind);
                ndat.nc().insert_add(tt);
                }
            }
        }
    dat.swap(ndat);
    is_.swap(nis_);
    return *this;
    }

IQTensor& IQTensor::
trace(const IQIndex& i1, const IQIndex& i2,
      const IQIndex& i3, const IQIndex& i4,
      const IQIndex& i5, const IQIndex& i6,
      const IQIndex& i7, const IQIndex& i8)
    {
    boost::array<IQIndex,NMAX> inds = {{ i1, i2, i3, i4,
                                i5, i6, i7, i8 }};
    trace(inds);
    return *this;
    }

/*
int IQTensor::
vecSize() const
	{
    if(this->isNull()) return 0;
	int s = 0;
    Foreach(const ITensor& t, dat())
	    s += t.vecSize();
	return s;
	}
    */

/*
//Redundant: use indices().dim() instead
int IQTensor::
maxSize() const
	{
    int ms = 1;
	for(int j = 0; j < is_->r(); ++j)
	    ms *= is_->m(j);
    return ms;
    }
*/

/*
void IQTensor::
assignToVec(VectorRef v) const
	{
	if(vecSize() != v.Length())
	    Error("Mismatched sizes in IQTensor::assignToVec(VectorRef v).");
	int off = 1;
    Foreach(const ITensor& t, blocks())
	    {
	    int d = t.vecSize();
	    t.assignToVec(v.SubVector(off,off+d-1));
	    off += d;
	    }
	}
    */

void IQTensor::
randomize() 
	{ 
    if(isNull())
        Error("Can't randomize default constructed IQTensor.");
    if(dat().empty())
        Error("Can't randomize IQTensor having no blocks");

	dat.solo(); 

    const QN D = div(*this);

    QCounter C(*is_);

    for(;C.notDone();++C)
        {
        QN nd;
        for(int n = 1; n <= r(); ++n)
            {
            nd += is_->index(n).dir()*is_->index(n).qn(1+C.i[n]);
            }
        if(nd != D) continue;

        IndexSet<Index> nset;
        for(int n = 1; n <= r(); ++n)
            {
            nset.addindex(is_->index(n).index(1+C.i[n]));
            }

        if(dat().hasBlock(nset))
            {
            dat.nc().get(nset).randomize();
            }
        else
            {
            ITensor t(nset);
            t.randomize();
            dat.nc().insert_add(t);
            }
        }
	}


IQTensor& IQTensor::
conj()
    {
    if(!this->isComplex())
        {
        soloIndex();
        is_->conj();
        return *this;
        }
    else
        {
        solo();
        is_->conj();
        Foreach(ITensor& t, dat.nc())
            {
            t.conj();
            }
        }
    return *this;
    }

void IQTensor::
swap(IQTensor& other)
    {
    is_.swap(other.is_);
    dat.swap(other.dat);
    }

std::ostream& 
operator<<(std::ostream & s, const IQTensor& T)
    {
	s << "/--------------IQTensor--------------\n";
    if(T.isNull())
        {
        s << "     (IQTensor is null)\n\n";
        return s;
        }
    s << "IQIndices:\n" << T.indices();
    s << "ITensor Blocks:\n";
    Foreach(const ITensor& t, T.blocks())
        { 
        if(t.r() > 0)
            {
            s << "  ";
            //Treat first Index specially in order to add trailing commas
            IndexSet<Index>::const_iterator it = t.indices().begin();
            const IQIndex& I1 = findIQInd(T,*it);
            s << it->name() << ":" << qn(I1,*it) << "<" << I1.dir() << ">";
            for(++it; it != t.indices().end(); ++it)
                {
                s << ", ";
                const IQIndex& I = findIQInd(T,*it);
                s << it->name() << ":" << qn(I,*it) << "<" << I.dir() << ">";
                }
            s << endl;
            }
        s << "  " << t << endl; 
        }
	s << "\\------------------------------------\n\n";
    return s;
    }

//Helper method for operator*= and operator/=
bool static
vectorARContains(const vector<ApproxReal>& v, ApproxReal r)
    {
    Foreach(const ApproxReal& x, v)
        {
        if(r == x) return true;
        }
    return false;
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

    solo();

    vector<ApproxReal> common_inds;
    
    //Load iqindex_ with those IQIndex's *not* common to *this and other
    boost::array<IQIndex,NMAX> riqind_holder;
    int rholder = 0;

    typedef IndexSet<IQIndex>::const_iterator
    const_iqind_it;

    for(int i = 1; i <= is_->r(); ++i)
        {
        const IQIndex& I = is_->index(i);
        const_iqind_it f = find(other.is_->begin(),other.is_->end(),I);
        if(f != other.is_->end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Global::checkArrows())
                {
                if(f->dir() == I.dir())
                    {
                    Print(this->indices());
                    Print(other.indices());
                    cout << "IQIndex from *this = " << I << endl;
                    cout << "IQIndex from other = " << *f << endl;
                    cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                    throw ArrowError("Incompatible arrow directions in IQTensor::operator*=.");
                    }
                }

            Foreach(const Index& i, I.indices())
                { 
                common_inds.push_back(i.uniqueReal()); 
                }

            common_inds.push_back(I.uniqueReal());
            }
        else 
            { 
            riqind_holder[rholder] = I;
            ++rholder;
            }
        }

    for(int i = 1; i <= other.is_->r(); ++i)
        {
        const IQIndex& I = other.is_->index(i);
        if(!vectorARContains(common_inds,I.uniqueReal()))
            { 
            if(rholder >= NMAX)
                {
                Print(this->indices());
                Print(other.indices());
                cout << "Uncontracted IQIndices found so far:" << endl;
                for(int n = 0; n < rholder; ++n)
                    {
                    cout << riqind_holder[n] << endl;
                    }
                Error("Too many indices (>= 8) on resulting IQTensor");
                }
            riqind_holder[rholder] = I;
            ++rholder;
            }
        }

    is_ = boost::make_shared<IndexSet<IQIndex> >(riqind_holder,rholder,0);

    IQTDat<ITensor>::StorageT old_itensor; 
    dat.nc().swap(old_itensor);

    typedef IQTDat<ITensor>::const_iterator
    cbit;

    typedef pair<ApproxReal,cbit>
    blockpair;

    vector<blockpair> other_block;
    other_block.reserve(other.dat().size());

    for(cbit ot = other.dat().begin(); ot != other.dat().end(); ++ot)
        {
        Real r = 0.0;
        Foreach(const Index& I, ot->indices())
            {
            if(vectorARContains(common_inds,I.uniqueReal()))
                r += I.uniqueReal(); 
            }
        other_block.push_back(make_pair(ApproxReal(r),ot));
        }

    ITensor prod;

    Foreach(const ITensor& t, old_itensor)
        {
        ApproxReal r(0);
        Foreach(const Index& I, t.indices())
            {
            if(vectorARContains(common_inds,I.uniqueReal()))
                r += I.uniqueReal();
            }
        Foreach(const blockpair& p, other_block)
            {
            if(r != p.first) continue;
            prod = t;
            prod *= *(p.second);
            if(prod.scale().sign() != 0)
                dat.nc().insert_add(prod);
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

    vector<ApproxReal> common_inds;
    
    boost::array<IQIndex,NMAX> riqind_holder;
    int rholder = 0;

    typedef IndexSet<IQIndex>::const_iterator
    const_iqind_it;

    for(int i = 1; i <= is_->r(); ++i)
        {
        const IQIndex& I = is_->index(i);
        const_iqind_it f = find(other.is_->begin(),other.is_->end(),I);
        if(f != other.is_->end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Global::checkArrows())
                if(f->dir() != I.dir())
                    {
                    Print(this->indices());
                    Print(other.indices());
                    cout << "IQIndex from *this = " << I << endl;
                    cout << "IQIndex from other = " << *f << endl;
                    cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                    throw ArrowError("Incompatible arrow directions in IQTensor::operator/=.");
                    }

            Foreach(const Index& i, I.indices())
                { 
                common_inds.push_back(i.uniqueReal()); 
                }

            common_inds.push_back(I.uniqueReal());
            }
        riqind_holder[rholder] = I;
        ++rholder;
        }

    bool inds_from_other = false;
    for(int i = 1; i <= other.is_->r(); ++i)
        {
        const IQIndex& I = other.is_->index(i);
        if(!vectorARContains(common_inds,I.uniqueReal()))
            { 
            if(rholder >= NMAX)
                {
                Print(this->indices());
                Print(other.indices());
                cout << "Uncontracted IQIndices found so far:" << endl;
                for(int n = 0; n < rholder; ++n)
                    {
                    cout << riqind_holder[n] << endl;
                    }
                Error("Too many indices (>= 8) on resulting IQTensor");
                }
            riqind_holder[rholder] = I;
            ++rholder;
            inds_from_other = true;
            }
        }

    //Only update IQIndices if they are different
    //from current set
    if(inds_from_other)
        {
        is_ = boost::make_shared<IndexSet<IQIndex> >(riqind_holder,rholder,0);
        }

    dat.solo();

    IQTDat<ITensor>::StorageT old_itensor; 
    dat.nc().swap(old_itensor);

    typedef IQTDat<ITensor>::const_iterator
    cbit;

    typedef pair<ApproxReal,cbit>
    blockpair;

    vector<blockpair> other_block;
    other_block.reserve(other.dat().size());

    for(cbit ot = other.dat().begin(); ot != other.dat().end(); ++ot)
        {
        Real r = 0.0;
        Foreach(const Index& I, ot->indices())
            {
            if(vectorARContains(common_inds,I.uniqueReal()))
                r += I.uniqueReal(); 
            }
        other_block.push_back(make_pair(ApproxReal(r),ot));
        }

    ITensor prod;

    Foreach(const ITensor& t, old_itensor)
        {
        ApproxReal r(0);
        Foreach(const Index& I, t.indices())
            {
            if(vectorARContains(common_inds,I.uniqueReal()))
                r += I.uniqueReal();
            }
        Foreach(const blockpair& p, other_block)
            {
            if(r != p.first) continue;
            prod = t;
            prod /= *(p.second);
            if(prod.scale().sign() != 0)
                dat.nc().insert_add(prod);
            }
        }

    return *this;

    } //IQTensor& IQTensor::operator/=(const IQTensor& other)

//Extracts the real and imaginary parts of the 
//component of a rank 0 tensor (scalar)
Complex IQTensor::
toComplex() const
    {
    if(this->isComplex())
        {
        Real re = realPart(*this).toReal();
        Real im = imagPart(*this).toReal();
        return Complex(re,im);
        }
    return Complex(toReal(),0);
    }

Real IQTensor::
toReal() const
    {
    if(is_->r() != 0)
        Error("IQTensor not a real scalar");
#ifdef DEBUG
    if(blocks().size() > 1)
        Error("Too many blocks");
#endif
    if(empty())
        return 0;
    else
        return dat().begin()->toReal();
    }

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    //TODO: account for fermion sign here

    if(this->isNull())
        {
        operator=(other);
        return *this;
        }

    if(this == &other) 
        {
        operator*=(2);
        return *this;
        }

    /*
    //EMS Mar 7 2013: not sure what this does or if it's correct
    if(is_->r() == 0)	// Automatic initializing a summed IQTensor in a loop
        { 
        return (*this = other); 
        }
        */

    IQTensor& This = *this;

    if(fabs(This.is_->uniqueReal()-other.is_->uniqueReal()) > 1.0e-11) 
        {
        Print(This.indices());
        Print(other.indices());
        Print(This.is_->uniqueReal());
        Print(other.is_->uniqueReal());
        Error("Mismatched indices in IQTensor::operator+=");
        }

    dat.solo(); 

    Foreach(const ITensor& t, other.dat())
        { 
        dat.nc().insert_add(t);
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
    if(this->isNull()) return ITensor();

    //if(Global::debug1())
    //    {
    //    cout << "Converting an IQTensor to ITensor" << endl;
    //    PAUSE
    //    }

    //Resulting ITensor's indices are 
    //the Index versions of this's IQIndices
    IndexSet<Index> indices;
    for(int j = 1; j <= is_->r(); ++j)
        {
        indices.addindex(Index(is_->index(j)));
        }
    ITensor res(indices);

    //Loop over ITensors (blocks) within this IQTensor
    Foreach(const ITensor& t, dat())
        {
        ITensor exp(t);
        //Loop over Index's of the k'th ITensor
        Foreach(const Index& small, t.indices())
            {
            //Want to transform 'small' into the 
            //Index version of the IQIndex that contains
            //it, with the appropriate offset

            //Find the IQIndex that contains 'small'
            const IQIndex* big = 0;
            Foreach(const IQIndex& I, this->indices())
                if(hasindex(I,small))
                    {
                    big = &I;
                    break;
                    }
            exp.expandIndex(small,*big,offset(*big,small));
            }
        //Once all Indices expanded, add to res
        res += exp;
        }
    return res;
    } //IQTensor::operator ITensor() const

void IQTensor::
soloIndex()
	{
	if(!is_)
        Error("IQTensor is null");

	if(!is_.unique())
        is_ = boost::make_shared<IndexSet<IQIndex> >(*is_);
    }


void IQTensor::
solo()
    {
    soloIndex();
    dat.solo();
    }

IQTensor& IQTensor::
takeRealPart()
    {
    dat.solo();
    Foreach(ITensor& t, dat.nc())
        {
        t.takeRealPart();
        }
    return *this;
    }

IQTensor& IQTensor::
takeImagPart()
    {
    dat.solo();
    Foreach(ITensor& t, dat.nc())
        {
        t.takeImagPart();
        }
    return *this;
    }

void IQTensor::
pseudoInvert(Real cutoff)
    {
    dat.solo();
    Foreach(ITensor& t, dat.nc())
        {
        t.pseudoInvert();
        }
    }

void IQTensor::
replaceIndex(const IQIndex& oind,
             const IQIndex& nind,
             const OptSet& opts)
    { 
    if(opts.getBool("CheckArrows",true))
        {
        if(nind.dir() != dir(*this,oind))
            {
            Error("replaceIndex: arrow dir's don't match.");
            }
        }
    if(nind.nindex() != oind.nindex())
        {
        Error("replaceIndex: different number of blocks.");
        }
    solo(); 
    is_->replaceIndex(oind,nind); 
    Foreach(ITensor& t, dat.nc())
        {
        Foreach(const Index& i, t.indices())
            {
            const int j = findindex(oind,i);
            if(j != 0)
                {
                t.replaceIndex(i,nind.index(j));
                }
            }
        }
    }


Real 
Dot(IQTensor x, const IQTensor& y)
    {
    IQIndex I = commonIndex(x,y);
    if(I.dir() == dir(y.indices(),I))
        {
        x.conj();
        }
    x *= y;
    return x.toReal();
    }

Complex 
BraKet(IQTensor x, const IQTensor& y)
    {
    x.conj();
    x *= y;
    return x.toComplex();
    }


QN
div(const IQTensor& T, const OptSet& opts)
	{
	if(T.empty())
	    {   
        Print(T);
	    Error("IQTensor has no blocks");
	    }

    //Calculate divergence of first block
    QN div_;
    IQTDat<ITensor>::const_iterator it = T.blocks().begin();
    Foreach(const Index& i, it->indices())
        {
        div_ += qn(T,i)*dir(T,i);
        }

    bool fast = opts.getBool("Fast",false);
#ifdef DEBUG
    fast = false;
#endif
    if(fast) return div_;

    //Check that remaining blocks have same divergence
    for(++it; it != T.blocks().end(); ++it)
        {
        QN q;
        Foreach(const Index& i, it->indices())
            {
            q += qn(T,i)*dir(T,i);
            }
        if(q != div_)
            {
            Global::printdat() = true;
            cout << "\n-------------------------\n" << endl;
            cout << "div = " << div_ << "\n" << endl;
            cout << "div of this block = " << q << "\n" << endl;
            cout << "Offending ITensor = \n" << *it << "\n" << endl;
            Print(T.indices());
            cout << "\n-------------------------\n" << endl;
            it = T.blocks().begin();
            cout << "First ITensor = \n" << *it << "\n" << endl;
            Error("Inconsistent divergence of IQTensor block");
            }
        }

	return div_;
	}

const IQIndex&
findIQInd(const IQTensor& T, const Index& i)
    {
    Foreach(const IQIndex& J, T.indices())
        {
        if(hasindex(J,i)) 
            return J;
        }
    Print(T.indices());
    Print(i);
    Error("Index i not found in any of T's IQIndices");
    return IQIndex::Null();
    }

QN
qn(const IQTensor& T, const Index& i)
	{
	return qn(findIQInd(T,i),i);
	} 

Arrow
dir(const IQTensor& T, const Index& i)
	{
    return findIQInd(T,i).dir();
	}

Arrow
dir(const IQTensor& T, const IQIndex& I)
	{
    Foreach(const IQIndex& J, T.indices())
        {
        if(I == J) return J.dir();
        }
    Error("dir: IQIndex not found");
    return Out;
	}

bool
uses_ind(const IQTensor& T, const Index& ii)
    {
    Foreach(const ITensor& t, T.blocks())
        {
        if(hasindex(t,ii)) 
            return true;
        }
    return false;
    }

bool
isZero(const IQTensor& T, const OptSet& opts)
    {
    if(T.empty()) return true;
    //done with all fast checks
    if(opts.getBool("Fast",false)) return false;
    Foreach(const ITensor& t, T.blocks())
        {
        if(!isZero(t)) return false;
        }
    return true;
    }
