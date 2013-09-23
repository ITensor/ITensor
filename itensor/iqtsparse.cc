//
// Distributed under the ITensor Library License, Version 1.0.
//    (See accompanying LICENSE file.)
//
#include "iqtsparse.h"
using namespace std;
using boost::format;
using boost::array;

//
// IQTSparse
//

IQTSparse::
IQTSparse()
    :
    is_(IndexSet<IQIndex>::Null()),
    d_(IQTDat<ITSparse>::Null())
    { }

IQTSparse::
IQTSparse(const IQIndex& i1)
    :
    is_(boost::make_shared<IndexSet<IQIndex> >(i1)),
    d_(boost::make_shared<IQTDat<ITSparse> >())
    { 
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2)
    :
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2)),
    d_(boost::make_shared<IQTDat<ITSparse> >())
    { 
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2,
          Real r)
    :
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2)),
    d_(boost::make_shared<IQTDat<ITSparse> >())
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
IQTSparse(const IQIndex& i1, const IQIndex& i2,
          const VectorRef& D)
    :
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2)),
    d_(boost::make_shared<IQTDat<ITSparse> >())
    { 
#ifdef DEBUG
    if(i1.m() != i2.m() || i1.nindex() != i2.nindex())
        {
        Print(i1);
        Print(i2);
        Error("IQIndices must have same size and number of blocks");
        }
    if(D.Length() != i1.m())
        {
        Error("Incorrect size of Vector");
        }
#endif
    int n = 1;
    for(int j = 1; j <= i1.nindex(); ++j)
        {
        int bsize = i1.index(j).m();
        operator+=(ITSparse(i1.index(j),i2.index(j),D.SubVector(n,n-1+bsize)));
        n += bsize;
        }
    }

IQTSparse::
IQTSparse(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3)
    :
    is_(boost::make_shared<IndexSet<IQIndex> >(i1,i2,i3)),
    d_(boost::make_shared<IQTDat<ITSparse> >())
    { 
    }

bool IQTSparse::
isNull() const { return (d_ == IQTDat<ITSparse>::Null() || is_ == IndexSet<IQIndex>::Null()); }

bool IQTSparse::
isComplex() const
    {
    Foreach(const ITSparse& t, blocks())
        {
        if(t.isComplex()) return true;
        }
    return false;
    }


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

IQTSparse& IQTSparse::
noprime(IndexType type)
    { 
    solo();

    is_->noprime(type); 

    Foreach(ITSparse& t, ncblocks())
        {
        t.noprime(type);
        }
    return *this;
    }

IQTSparse& IQTSparse::
prime(IndexType type, int inc) 
    { 
    solo();

    is_->prime(type,inc);

    Foreach(ITSparse& t, ncblocks())
        {
        t.prime(type,inc);
        }
    return *this;
    }

IQTSparse& IQTSparse::
mapprime(int plevold, int plevnew, IndexType type)
    { 
    solo();

    is_->mapprime(plevold,plevnew,type); 

    Foreach(ITSparse& t, ncblocks())
        {
        t.mapprime(plevold,plevnew,type);
        }
    return *this;
    }

IQTSparse& IQTSparse::
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
    return *this;
    }

IQTSparse& IQTSparse::
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
    return *this;
    }

void IQTSparse::
conj()
    {
    if(!this->isComplex())
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
scaleOutNorm() 
	{
    Real nrm = norm();
    ncblocks().scaleTo(nrm);
	}

void IQTSparse::
scaleTo(const LogNumber& newscale)
	{
    ncblocks().scaleTo(newscale);
	}


void IQTSparse::
read(std::istream& s)
    {
    bool null_;
    s.read((char*) &null_,sizeof(null_));
    if(null_) 
        { *this = IQTSparse(); return; }
    is_ = boost::make_shared<IndexSet<IQIndex> >();
    is_->read(s);
    d_ = boost::make_shared<IQTDat<ITSparse> >();
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
        const IQTDat<ITSparse>& old_dat(*d_);
        d_ = boost::make_shared<IQTDat<ITSparse> >();
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
        is_ = boost::make_shared<IndexSet<IQIndex> >();
        *is_ = old_is;
        }
    }

void IQTSparse::
solo()
    {
    soloIndex();
    soloDat();
    }

bool static
vectoruRContains(const vector<Real>& v, Real r)
    {
    Foreach(const Real& x, v)
        {
        if(fabs(r - x) <= UniqueRealAccuracy) return true;
        }
    return false;
    }

void
product(const IQTSparse& S, const IQTensor& T, IQTensor& res)
    {
    if(S.isNull()) 
        Error("IQTSparse null in product");

    if(T.isNull()) 
        Error("Multiplying by null IQTensor");

    vector<Real> common_inds;
    
    //Load iqindex_ with those IQIndex's *not* common to *this and other
    static vector<IQIndex> riqind_holder;
    riqind_holder.resize(0);

    for(int i = 1; i <= S.is_->r(); ++i)
        {
        const IQIndex& I = S.is_->index(i);
        IndexSet<IQIndex>::const_iterator f = find(T.is_->begin(),T.is_->end(),I);
        if(f != T.is_->end()) //I is an element of other.iqindex_
            {
            //Check that arrow directions are compatible
            if(Global::checkArrows())
                if(f->dir() == I.dir())
                    {
                    Print(S.indices());
                    Print(T.indices());
                    cerr << "IQIndex from S = " << I << endl;
                    cerr << "IQIndex from T = " << *f << endl;
                    cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                    throw ArrowError("Incompatible arrow directions in IQTensor::operator*=.");
                    }

            Foreach(const Index& i, I.indices())
                { common_inds.push_back(i.uniqueReal()); }

            common_inds.push_back(I.uniqueReal());
            }
        else 
            { 
            riqind_holder.push_back(I); 
            }
        }

    for(int i = 1; i <= T.is_->r(); ++i)
        {
        const IQIndex& I = T.is_->index(i);
        if(!vectoruRContains(common_inds,I.uniqueReal()))
            { 
            riqind_holder.push_back(I); 
            }
        }

    res = IQTensor(riqind_holder);

    typedef IQTDat<ITSparse>::const_iterator
    cbit;

    typedef pair<Real,cbit>
    blockpair;

    vector<blockpair> Sblock;
    Sblock.reserve(S.blocks().size());

    for(cbit ot = S.blocks().begin(); ot != S.blocks().end(); ++ot)
        {
        Real r = 0.0;
        Foreach(const Index& I, ot->indices())
            {
            if(vectoruRContains(common_inds,I.uniqueReal()))
                r += I.uniqueReal(); 
            }
        Sblock.push_back(make_pair(r,ot));
        }

    ITensor prod;

    Foreach(const ITensor& t, T.blocks())
        {
        Real r = 0;
        Foreach(const Index& I, t.indices())
            {
            if(vectoruRContains(common_inds,I.uniqueReal()))
                r += I.uniqueReal();
            }
        Foreach(const blockpair& p, Sblock)
            {
            if(fabs(r - p.first) > UniqueRealAccuracy) continue;
            prod = t;
            prod *= *(p.second);
            if(prod.scale().sign() != 0)
                res.dat.nc().insert_add(prod);
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
