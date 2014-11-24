//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "iqtensor.h"
#include "qcounter.h"
#include <algorithm>

namespace itensor {

using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::find;
using std::pair;
using std::make_pair;

void static
insertAdd(ITensor& lhs, const ITensor& rhs)
    {
    if(!lhs.valid()) lhs = rhs;
    else             lhs += rhs;
    }

int static
indPos(const Index& i,
       const IQIndex& J)
    {
    for(int j = 0; j < J.nindex(); ++j)
        {
        if(J[j] == i) return j;
        }
    return -1;
    }

// Determine the position of a given block (specified by an IndexSet<Index>, 
// which could be from an ITensor or a set of IQIndexVals)
// in the storage of an IQTensor (specified by an IndexSet<IQIndex>)
int static
blockPos(const IndexSet<Index>& is,
         const IndexSet<IQIndex>& qs)
    {
    int pos = 0;
    int dim = 1;

    for(const IQIndex& J : qs)
        {
        bool found = false;
        for(const Index& i : is)
            {
            if(hasindex(J,i))
                {
                pos += dim*indPos(i,J);
                dim *= J.nindex();
                found = true;
                break;
                }
            }
        if(!found) throw ITError("Index not found in IndexSet<IQIndex>");
        }

    return pos;
    }


//
// IQTensor
//

bool IQTensor::
valid() const 
    { 
    return bool(d_) && (d_.get() != IQTDat::Null().get());
    }

bool IQTensor::
isComplex() const
    {
    for(const ITensor& t : *d_)
        {
        if(t.isComplex()) return true;
        }
    return false;
    }

int IQTensor::
r() const 
    { 
    return is_.r(); 
    }

bool IQTensor::
empty() const { return d_->empty(); }

//----------------------------------------------------
//IQTensor: Constructors 

void IQTensor::
allocate()
    {
    int dim = 1;
    for(const IQIndex& I : this->indices())
        {
        dim *= I.nindex();
        }
    d_ = make_shared<Storage>(dim);
    }


IQTensor::
IQTensor() 
    :
    d_(IQTDat::Null())
    { }

IQTensor::
IQTensor(Real val) 
    { 
    allocate();
    operator+=(ITensor(val));
    }

IQTensor::
IQTensor(const IQIndex& i1) 
    : 
    is_(i1),
    d_(make_shared<Storage>(i1.nindex()))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2) 
    : 
    is_(i1,i2),
    d_(make_shared<Storage>(i1.nindex()*i2.nindex()))
    { 
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) 
	: 
    is_(i1,i2,i3)
    { 
    allocate();
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4) 
    : 
    is_(i1,i2,i3,i4)
    { 
    allocate();
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5) 
    : 
    is_(i1,i2,i3,i4,i5)
    { 
    allocate();
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6)
	: 
    is_(i1,i2,i3,i4,i5,i6)
	{ 
    allocate();
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7)
	: 
    is_(i1,i2,i3,i4,i5,i6,i7)
	{ 
    allocate();
    }

IQTensor::
IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
         const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
         const IQIndex& i7,const IQIndex& i8)
	: 
    is_(i1,i2,i3,i4,i5,i6,i7,i8)
	{ 
    allocate();
    }

IQTensor::
IQTensor(vector<IQIndex>& iqinds_) 
	: 
    is_(iqinds_)
	{ 
#ifdef DEBUG
    for(const IQIndex& I : iqinds_)
        {
        if(I == IQIndex::Null())
            Error("IQIndex is null");
        }
#endif
    allocate();
    }

IQTensor::
IQTensor(const IQIndexVal& iv1) 
    : 
    is_(iv1.index),
    d_(make_shared<Storage>(iv1.index.nindex()))
	{ 
	operator()(iv1) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2) 
	: 
    is_(iv1.index,iv2.index),
    d_(make_shared<Storage>(iv1.index.nindex()*iv2.index.nindex()))
	{ 
    operator()(iv1,iv2) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
         const IQIndexVal& iv3)
	: 
    is_(iv1.index,iv2.index,iv3.index)
	{ 
    allocate();
    operator()(iv1,iv2,iv3) = 1;
	}

IQTensor::
IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
         const IQIndexVal& iv3, const IQIndexVal& iv4)
	: 
    is_(iv1.index,iv2.index,iv3.index,iv4.index)
	{ 
    allocate();
    operator()(iv1,iv2,iv3,iv4) = 1;
	}

void IQTensor::
read(std::istream& s)
    {
    bool null_;
    s.read((char*) &null_,sizeof(null_));
    if(null_) 
        { *this = IQTensor(); return; }
    is_.read(s);
    d_ = make_shared<Storage>();
    d_->read(s);
    }

void IQTensor::
write(std::ostream& s) const
	{
	bool null_ = !valid();
	s.write((char*) &null_,sizeof(null_));
	if(null_) return;
    is_.write(s);
	d_->write(s);
	}

IQTensor& IQTensor::
operator*=(Real fac) 
    { 
    solo();

    if(fac == 0) 
        { 
        d_->clear(); 
        return *this; 
        }

    for(ITensor& t : *d_)
        {
        t *= fac;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator/=(Real fac) 
    { 
    solo();

    for(ITensor& t : *d_)
        {
        t /= fac;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator*=(Complex z) 
    { 
    solo();

    if(z.real() == 0 && z.imag() == 0) 
        { 
        d_->clear();
        return *this; 
        }

    for(ITensor& t : *d_)
        {
        t *= z;
        }

    return *this; 
    }

IQTensor& IQTensor::
operator*=(const LogNumber& lgnum) 
    { 
    solo();

    for(ITensor& t : *d_)
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
        solo();
        getBlock(t.indices()) = t;
        }
    }

IQTensor& IQTensor::
operator+=(const ITensor& t) 
    { 
#ifdef DEBUG
    if(!this->valid())
        Error("Calling operator+=(ITensor) on an invalid IQTensor.");
    if(!this->empty())
        {
        const
        QN d = div(*this);
        QN q;
        for(const Index& i : t.indices())
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
        solo();
        insertAdd(getBlock(t.indices()),t);
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
    solo();
    array<IQIndexVal,NMAX+1> iv 
        = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};

    int nn = 0; 
    IndexSet<Index> is;
    while(GET(iv,nn+1) != IQIndexVal::Null()) 
        {
        ++nn;
        if(!hasindex(*this,iv.at(nn).index)) 
            Error("IQTensor::operator(): IQIndex not found.");
        is.addindex(iv[nn].indexqn());
        }
    if(nn != r()) 
        {
        Error("Wrong number of IQIndexVals provided");
        }

    ITensor& block = getBlock(is);
    if(!block.valid()) block = ITensor(is);

    return block.operator()(iv1.blockIndexVal(),
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
    array<IQIndexVal,NMAX+1> iv 
        = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};

    int nn = 0; 
    IndexSet<Index> is;
    while(GET(iv,nn+1) != IQIndexVal::Null()) 
        {
        //ur += GET(iv,++nn).indexqn().uniqueReal(); 
        if(!hasindex(*this,iv.at(nn).index)) 
            Error("IQTensor::operator(): IQIndex not found.");
        is.addindex(iv[nn].indexqn());
        }
    if(nn != r()) 
        Error("Wrong number of IQIndexVals provided");

    const ITensor& block = getBlock(is);
    if(!block)
        {
        return 0.;
        }
    else
        {
        return block.operator()(iv1.blockIndexVal(),
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

    is_.noprime(type);

    for(ITensor& t : *d_)
        t.noprime(type); 
    return *this;
	} 

IQTensor& IQTensor::
prime(IndexType type, int inc)
	{
	solo();

    is_.prime(type,inc);

    for(ITensor& t : *d_)
	    t.prime(type,inc);
    return *this;
	}

IQTensor& IQTensor::
mapprime(int plevold, int plevnew, IndexType type)
    {
    solo();

    is_.mapprime(plevold,plevnew,type);

    for(ITensor& t : *d_)
	    t.mapprime(plevold,plevnew,type);
    return *this;
	}

IQTensor& IQTensor::
prime(const IQIndex& I, int inc)
	{
	solo();

    is_.prime(I,inc);

    for(ITensor& t : *d_)
    for(const Index& i : I.indices())
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

    is_.noprime(I);

    for(ITensor& t : *d_)
    for(const Index& i : I.indices())
        {
        if(hasindex(t,i)) 
            t.noprime(i);
        }
    return *this;
	}




LogNumber IQTensor::
normLogNum() const
    {
    if(d_->empty()) return LogNumber(0);

    Real maxLogNum = -maxlogdouble;
    for(const ITensor& t : *d_)
        { 
        if(t.scale().logNum() > maxLogNum)
            maxLogNum = t.scale().logNum();
        }

    //Add sum of squares of the block norms with exp(2*maxLogNum) scaled out, 
    //using lognorm as a temporary
    Real lognorm = 0;
    for(const ITensor& t : *d_)
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

Vector IQTensor::
diag() const
    {
    if(this->isComplex()) 
        Error("diag() may only be called on real IQTensors - try taking real or imaginary part first");
    int nb = this->indices().front().nindex();
    for(const IQIndex& I : this->indices())
        nb = min(nb,I.nindex());
    vector<Real> els;
    for(int n = 1; n <= nb; ++n)
        {
        IndexSet<Index> is;
        int bsize = this->indices().front().index(n).m();
        for(const IQIndex& I : this->indices())
            {
            is.addindex(I.index(n));
            bsize = min(bsize,I.index(n).m());
            }

        Vector d;

        const ITensor& block = getBlock(is);
        if(block) d = block.diag();
        else      d = Vector(bsize,0);

        for(int j = 1; j <= d.Length(); ++j)
            {
            els.push_back(d(j));
            }
        }
    Vector D(els.size());
    for(int j = 0; j < D.Length(); ++j)
        {
        D[j] = els[j];
        }
    return D;
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
    for(const ITensor& t : T.blocks())
        { 
        sum += sumels(t); 
        }
    return sum;
    }

void IQTensor::
scaleOutNorm()
    {
    solo(); 
    const LogNumber newscale = normLogNum();
    for(ITensor& t : *d_)
        t.scaleTo(newscale);
    }

void IQTensor::
scaleTo(const LogNumber& newscale)
    {
    solo(); 
    for(ITensor& t : *d_)
        t.scaleTo(newscale);
    }

void IQTensor::
clean(Real min_norm)
    { 
    solo(); 
    for(ITensor& t : *d_)
        {
        if(t.norm() < min_norm)
            t = ITensor();
        }
    }


void IQTensor::
tieIndices(const array<IQIndex,NMAX>& indices, 
           int niqind, 
           const IQIndex& tied)
    {
    if(niqind < 1) Error("No IQIndices to tie");

    const int nindex = indices[0].nindex();

    IndexSet<IQIndex> nis_(tied);

    int nmatched = 0;
    for(int k = 1; k <= is_.r(); ++k)
        {
        const IQIndex& K = is_.index(k);
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
            nis_.addindex(K);
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

    is_.swap(nis_);

    IQTDatPtr prevdat;
    prevdat.swap(d_);
    allocate();

    array<Index,NMAX> totie;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 0; n < niqind; ++n)
            totie[n] = indices[n].index(i);

        for(const ITensor& t : *prevdat)
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
                insertAdd(getBlock(nt.indices()),nt);
                }
            }
        }
    }

void IQTensor::
tieIndices(const IQIndex& i1, const IQIndex& i2, const IQIndex& tied)
    {
    array<IQIndex,NMAX> inds =
        {{ i1, i2, 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null(), 
           IQIndex::Null(), IQIndex::Null() }};

    tieIndices(inds,2,tied);
    }

IQTensor& IQTensor::
trace(const array<IQIndex,NMAX>& indices, int niqind)
    {
    if(niqind < 0)
        {
        niqind = 0;
        while(indices[niqind] != IQIndex::Null()) ++niqind;
        }

    if(niqind < 1) Error("No IQIndices to trace");

    const int nindex = indices[0].nindex();
    const int tm = indices[0].m();

    IndexSet<IQIndex> nis_;

    int nmatched = 0;
    for(int k = 1; k <= is_.r(); ++k)
        {
        const IQIndex& K = is_.index(k);
        bool K_traced = false;
        for(int j = 0; j < niqind; ++j)
            {
            if(K == indices[j]) 
                {
                if(indices[j].m() != tm)
                    Error("Traced indices must have matching m's");
                K_traced = true;
                ++nmatched;
                break;
                }
            }
        if(!K_traced)
            {
            nis_.addindex(K);
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

    is_.swap(nis_);

    IQTDatPtr prevdat;
    prevdat.swap(d_);
    allocate();

    array<Index,NMAX> totrace;
    for(int i = 1; i <= nindex; ++i)
        {
        for(int n = 0; n < niqind; ++n)
            {
            totrace[n] = indices[n].index(i);
            }

        for(const ITensor& t : *prevdat)
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
                insertAdd(getBlock(tt.indices()),tt);
                }
            }
        }


    return *this;
    }

IQTensor& IQTensor::
trace(const IQIndex& i1, const IQIndex& i2,
      const IQIndex& i3, const IQIndex& i4,
      const IQIndex& i5, const IQIndex& i6,
      const IQIndex& i7, const IQIndex& i8)
    {
    array<IQIndex,NMAX> inds = {{ i1, i2, i3, i4,
                                i5, i6, i7, i8 }};
    trace(inds);
    return *this;
    }

void IQTensor::
randomize(const Args& args) 
	{ 
    if(!valid())
        Error("Can't randomize default constructed IQTensor.");
    if(d_->empty())
        Error("Can't randomize IQTensor having no blocks");

	solo(); 

    const QN D = div(*this);

    QCounter C(is_);

    for(;C.notDone();++C)
        {
        QN nd;
        for(int n = 1; n <= r(); ++n)
            {
            nd += is_.index(n).dir()*is_.index(n).qn(1+C.i[n]);
            }
        if(nd != D) continue;

        IndexSet<Index> nset;
        for(int n = 1; n <= r(); ++n)
            {
            nset.addindex(is_.index(n).index(1+C.i[n]));
            }

        ITensor& block = getBlock(nset);
        if(block)
            {
            block.randomize(args);
            }
        else
            {
            ITensor t(nset);
            t.randomize(args);
            insertAdd(block,t);
            }
        }
	}

IQTensor& IQTensor::
conj()
    {
    if(!this->isComplex()) return *this;

    solo();
    for(ITensor& t : *d_)
        {
        t.conj();
        }
    return *this;
    }

IQTensor& IQTensor::
dag()
    {
    if(!this->isComplex())
        {
        is_.dag();
        return *this;
        }
    else
        {
        solo();
        is_.dag();
        for(ITensor& t : *d_)
            {
            t.dag();
            }
        }
    return *this;
    }

void IQTensor::
swap(IQTensor& other)
    {
    is_.swap(other.is_);
    d_.swap(other.d_);
    }

std::ostream& 
operator<<(std::ostream & s, const IQTensor& T)
    {
	s << "/--------------IQTensor--------------\n";
    if(!T)
        {
        s << "     (IQTensor is null)\n\n";
        return s;
        }
    s << "IQIndices:\n" << T.indices();
    s << "ITensor Blocks:\n";
    for(const ITensor& t : T.blocks())
        { 
        if(t.r() > 0)
            {
            s << "  ";
            //s << blockPos(t.indices(),T.indices()) << "  ";

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

int
dot_(const array<const int*,NMAX>& i,
     const array<int,NMAX>&  L)
    {
    int d = 0;
    array<const int*,NMAX>::const_iterator it = i.begin();
    array<int,NMAX>::const_iterator  Lt = L.begin();
    for(; Lt != L.end(); ++Lt,++it)
        {
        d += (*(*it)) * (*Lt);
        }
    return d;
    }


// Contracting product

struct BlockInfo
    {
    int u, //"partial index" associated with uncontracted indices
        c, //"partial index" associated with contracted indices
        p; //position of block in storage

    BlockInfo(int n) : u(0),c(0),p(n) { }
    };
std::ostream&
operator<<(std::ostream& s, BlockInfo nfo)
    {
    return s << "(" << nfo.u << "," << nfo.c << "," << nfo.p << ")";
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

    if(!this->valid()) 
        Error("'This' IQTensor null in product");

    if(!other)
        Error("Multiplying by null IQTensor");

    array<bool,NMAX> contractedL,
                     contractedR;
    contractedL.fill(false);
    contractedR.fill(false);

    //cdL/R is "contracted dimensions": weights/shapes 
    //of indices enumerating the contracted indices
    array<int,NMAX> cdL,
                    cdR;
    cdL.fill(-1000);
    cdR.fill(-1000);

    for(int i = 0, cdim = 1; i < is_.rn(); ++i)
        {
        const IQIndex& I = is_[i];
        for(int j = 0; j < other.is_.rn(); ++j)
            {
            const IQIndex& J = other.is_[j];
            if(J == I)
                {
                //I(==J) is contracted:

                //Check that arrow directions are compatible
                if(Global::checkArrows())
                    {
                    if(J.dir() == I.dir())
                        {
                        Print(this->indices());
                        Print(other.indices());
                        cout << "IQIndex from *this = " << I << endl;
                        cout << "IQIndex from other = " << J << endl;
                        cout << "Incompatible arrow directions in IQTensor::operator*=" << endl;
                        throw ArrowError("Incompatible arrow directions in IQTensor::operator*=.");
                        }
                    }

                contractedL[i] = true;
                contractedR[j] = true;

                cdL[i] = cdim;
                cdR[j] = cdim;
                cdim *= J.nindex();

                break;
                }
            }
        }

    //Finish making contractedL/R for m==1 IQIndices
    for(int i = is_.rn(); i < is_.r(); ++i)
    for(int j = other.is_.rn(); j < other.is_.r(); ++j)
        {
        if(other.is_[j] == is_[i]) 
            {
            contractedL[i] = true;
            contractedR[j] = true;
            break;
            }
        }

    //Load newindex with those IQIndex's *not* common to *this and other
    array<IQIndex,NMAX> newindex;
    int nnew = 0; //number of indices of product
    int nsize = 1; //size of storage for product

    //ud is "uncontracted dimensions": weights/shapes 
    //of indices enumerating the uncontracted indices
    array<int,NMAX> ud;
    ud.fill(-1000);
    int udim = 1; //accumulates products of dims of uncontracted indices

    for(int i = 0; i < is_.rn(); ++i)
        if(!contractedL[i])
            {
            const IQIndex& I = is_[i];
            newindex[nnew] = I;
            nsize *= I.nindex();
            ud[nnew] = udim;
            udim *= I.nindex();
            ++nnew;
            }

    //nucL is number of uncontracted from Left tensor
    const int nucL = nnew;

    for(int j = 0; j < other.is_.rn(); ++j)
        if(!contractedR[j])
            {
            const IQIndex& J = other.is_[j];
            newindex[nnew] = J;
            nsize *= J.nindex();
            ud[nnew] = udim;
            udim *= J.nindex();
            ++nnew;
            }

    //Finish making newindex by appending uncontracted m==1 IQIndices:
    for(int i = is_.rn(); i < is_.r(); ++i)
        if(!contractedL[i])
            {
            newindex[nnew++] = is_[i];
            }
    
    for(int j = other.is_.rn(); j < other.is_.r(); ++j)
        if(!contractedR[j])
            {
            newindex[nnew++] = other.is_[j];
            }


    IndexSet<IQIndex> nis(newindex,nnew,0);
    IQTDatPtr ndat = make_shared<IQTDat>(nsize);

    const IQTDat::Storage& L = d_->store();
    const IQTDat::Storage& R = other.d_->store();
    IQTDat::Storage& N = ndat->store();

    vector<BlockInfo> Lb;

    for(size_t p = 0; p < L.size(); ++p)
        {
        if(!L[p].valid()) continue;

        BlockInfo l(p);
        for(int n = 0, j = p, nu = 0; n < is_.rn(); ++n)
            {
            const int N = is_[n].nindex();
            const int i = j%N;
            j /= N;
            if(contractedL[n]) l.c += i*cdL[n];
            else               l.u += i*ud[nu++];
            }
        Lb.push_back(l);
        }

    for(size_t p = 0; p < R.size(); ++p)
        {
        if(!R[p].valid()) continue;

        BlockInfo r(p);
        for(int n = 0, j = p, nu = nucL; n < other.is_.rn(); ++n)
            {
            const int N = other.is_[n].nindex();
            const int i = j%N;
            j /= N;
            if(contractedR[n]) r.c += i*cdR[n];
            else               r.u += i*ud[nu++];
            }

        for(const BlockInfo& l : Lb)
            {
            if(l.c == r.c) insertAdd(N[l.u+r.u],L[l.p] * R[r.p]);
            }
        }

    is_.swap(nis);
    d_.swap(ndat);
    
    return *this;

    } //IQTensor& IQTensor::operator*=(const IQTensor& other)

IQTensor& IQTensor::
operator/=(const IQTensor& other)
    {
    //TODO: account for fermion sign here
    if(this == &other)
        {
        IQTensor cp_oth(other);
        return operator*=(cp_oth);
        }

    if(!this->valid()) 
        Error("'This' IQTensor null in product");

    if(!other)
        Error("Multiplying by null IQTensor");

    Counter u;

    const int zero = 0;

    array<const int*,NMAX> li,ri;
    for(size_t n = 0; n < li.size(); ++n)
        {
        li[n] = &zero;
        ri[n] = &zero;
        }

    //Load newindex with those IQIndex's *not* common to *this and other
    array<IQIndex,NMAX> newindex;
    int nnew = 0; //number of indices on product
    int nsize = 1; //size of storage for product

    for(int i = 0; i < is_.r(); ++i)
        {
        const IQIndex& I = is_[i];
        int j = 0;
        for(; j < other.is_.r(); ++j)
            {
            const IQIndex& J = other.is_[j];
            if(J == I)
                {
                //Check that arrow directions are compatible
                if(Global::checkArrows())
                    {
                    if(J.dir() != I.dir())
                        {
                        Print(this->indices());
                        Print(other.indices());
                        cout << "IQIndex from *this = " << I << endl;
                        cout << "IQIndex from other = " << J << endl;
                        cout << "Incompatible arrow directions in IQTensor::operator/=" << endl;
                        throw ArrowError("Incompatible arrow directions in IQTensor::operator/=");
                        }
                    }

                ++u.rn;
                u.n[u.rn] = J.nindex();
                li[i] = &(u.i[u.rn]);
                ri[j] = &(u.i[u.rn]);

                newindex[nnew] = J;
                nsize *= J.nindex();
                ++nnew;

                break;
                }
            }

        if(j == other.is_.r()) 
            { 
            // I is not contracted 
            ++u.rn;
            u.n[u.rn] = I.nindex();
            li[i] = &(u.i[u.rn]);
            newindex[nnew] = I;
            nsize *= I.nindex();
            ++nnew;
            }
        }

    for(int j = 0; j < other.is_.r(); ++j)
        {
        const IQIndex& J = other.is_[j];
        bool contracted = false;
        for(const IQIndex& I : this->indices())
            {
            if(I == J)
                {
                contracted = true;
                break;
                }
            }
        if(!contracted)
            {
            ++u.rn;
            u.n[u.rn] = J.nindex();
            ri[j] = &(u.i[u.rn]);
            newindex[nnew] = J;
            nsize *= J.nindex();
            ++nnew;
            }
        }

    array<int,NMAX> ll,
                    rl;
    std::fill(ll.begin(),ll.end(),0);
    std::fill(rl.begin(),rl.end(),0);

    for(int n = 0, dim = 1; n < is_.r(); ++n)
        {
        ll[n] = dim;
        dim *= is_[n].nindex();
        }
    for(int n = 0, dim = 1; n < other.is_.r(); ++n)
        {
        rl[n] = dim;
        dim *= other.is_[n].nindex();
        }

    is_ = IndexSet<IQIndex>(newindex,nnew,0);

    IQTDatPtr nd_ = make_shared<IQTDat>(nsize);

    const ITensor* pL = &(d_->store().front());
    const ITensor* pR = &(other.d_->store().front());
    ITensor* pN = &(nd_->store().front());

    for(; u.notDone(); ++u)
        {
        ITensor& n = pN[u.ind];
        const ITensor& l = pL[dot_(li,ll)];
        const ITensor& r = pR[dot_(ri,rl)];
        if(l.valid() && r.valid())
            {
            n = l;
            n /= r;
            }
        }

    d_.swap(nd_);

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
    if(is_.r() != 0)
        Error("IQTensor not a real scalar");
#ifdef DEBUG
    if(blocks().size() > 1)
        Error("Too many blocks");
#endif
    if(empty())
        return 0;
    else
        return d_->begin()->toReal();
    }

IQTensor& IQTensor::
operator+=(const IQTensor& other)
    {
    if(!this->valid())
        {
        Error("Calling += on null IQTensor");
        return *this;
        }
    if(!other.valid())
        {
        Error("Right-hand side of IQTensor += is null");
        return *this;
        }

    if(this == &other) 
        {
        operator*=(2);
        return *this;
        }

    solo(); 

    for(const ITensor& t : *(other.d_))
        { 
        insertAdd(getBlock(t.indices()),t);
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
    if(!this->valid()) return ITensor();

    //if(Global::debug1())
    //    {
    //    cout << "Converting an IQTensor to ITensor" << endl;
    //    PAUSE
    //    }

    //Resulting ITensor's indices are 
    //the Index versions of this's IQIndices
    IndexSet<Index> indices;
    for(int j = 1; j <= is_.r(); ++j)
        {
        indices.addindex(Index(is_.index(j)));
        }
    ITensor res(indices);

    //Loop over ITensors (blocks) within this IQTensor
    for(const ITensor& t : *d_)
        {
        ITensor exp(t);
        //Loop over Index's of the k'th ITensor
        for(const Index& small : t.indices())
            {
            //Want to transform 'small' into the 
            //Index version of the IQIndex that contains
            //it, with the appropriate offset

            //Find the IQIndex that contains 'small'
            const IQIndex* big = 0;
            for(const IQIndex& I : this->indices())
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


IQTensor& IQTensor::
takeRealPart()
    {
    if(isComplex())
        {
        solo();
        for(ITensor& t : *d_)
            {
            t.takeRealPart();
            }
        }
    return *this;
    }

IQTensor& IQTensor::
takeImagPart()
    {
    solo();
    for(ITensor& t : *d_)
        {
        t.takeImagPart();
        }
    return *this;
    }

void IQTensor::
pseudoInvert(Real cutoff)
    {
    solo();
    for(ITensor& t : *d_)
        {
        t.pseudoInvert();
        }
    }

void IQTensor::
replaceIndex(const IQIndex& oind,
             const IQIndex& nind,
             const Args& args)
    { 
    if(args.getBool("CheckArrows",true))
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
    is_.replaceIndex(oind,nind); 
    for(ITensor& t : *d_)
        {
        for(const Index& i : t.indices())
            {
            const int j = findindex(oind,i);
            if(j != 0)
                {
                t.replaceIndex(i,nind.index(j));
                }
            }
        }
    }

const ITensor& IQTensor::
getBlock(const IndexSet<Index>& inds) const
    {
    if(!valid()) Error("Default initialized IQTensor");
#ifdef DEBUG
    const int pos = blockPos(inds,is_);
    if(pos < 0 || pos >= d_->size()) 
        {
        Print(inds);
        Print(is_);
        Print(blockPos(inds,is_));
        Print(d_->size());
        Error("blockPos out of range");
        }
    const ITensor& t = d_->at(pos);
    return t;
#else
    return d_->at(blockPos(inds,is_));
#endif
    }

ITensor& IQTensor::
getBlock(const IndexSet<Index>& inds)
    {
#ifdef DEBUG
    if(!valid()) Error("Default initialized IQTensor");
#endif
    ITensor& t = d_->at(blockPos(inds,is_));
    return t;
    }

void IQTensor::
solo()
    {
#ifdef DEBUG
	if(!valid()) Error("IQTensor is default initialized");
#endif
	if(!d_.unique())
        {
        d_ = make_shared<IQTDat>(*d_);
        }
	}


Real 
Dot(IQTensor x, const IQTensor& y)
    {
    IQIndex I = commonIndex(x,y);
    if(I.dir() == dir(y.indices(),I))
        {
        x.dag();
        }
    x *= y;
    return x.toReal();
    }

Complex 
BraKet(IQTensor x, const IQTensor& y)
    {
    x.dag();
    x *= y;
    return x.toComplex();
    }


QN
div(const IQTensor& T, const Args& args)
	{
	if(T.empty())
	    {   
        Print(T);
	    Error("IQTensor has no blocks");
	    }

    //Calculate divergence of first block
    QN div_;
    IQTDat::const_iterator it = T.blocks().begin();
    for(const Index& i : it->indices())
        {
        div_ += qn(T,i)*dir(T,i);
        }

    bool fast = args.getBool("Fast",false);
#ifdef DEBUG
    fast = false;
#endif
    if(fast) return div_;

    //Check that remaining blocks have same divergence
    for(++it; it != T.blocks().end(); ++it)
        {
        QN q;
        for(const Index& i : it->indices())
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
    for(const IQIndex& J : T.indices())
        {
        if(hasindex(J,i)) 
            return J;
        }
    Print(T.indices());
    Print(i);
    throw ITError("Index i not found in any of T's IQIndices");
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
    for(const IQIndex& J : T.indices())
        {
        if(I == J) return J.dir();
        }
    Error("dir: IQIndex not found");
    return Out;
	}

bool
uses_ind(const IQTensor& T, const Index& ii)
    {
    for(const ITensor& t : T.blocks())
        {
        if(hasindex(t,ii)) 
            return true;
        }
    return false;
    }

bool
isZero(const IQTensor& T, const Args& args)
    {
    if(T.empty()) return true;
    //done with all fast checks
    if(args.getBool("Fast",false)) return false;
    for(const ITensor& t : T.blocks())
        {
        if(!isZero(t)) return false;
        }
    return true;
    }

void
checkStorage(const IQTensor& T)
    {
    const IQTDat& store = T.blocks();
    for(int n = 0; n < store.maxSize(); ++n)
        {
        const ITensor& b = store.at(n);
        if(b.valid())
            {
            const int pos = blockPos(b.indices(),T.indices());
            if(pos != n)
                {
                printfln("T.indices() = %s",T.indices());
                printfln("b = %s",b);
                printfln("pos=%d, n=%d",pos,n);
                throw ITError("ITensor block in wrong storage position.");
                }
            }
        }
    }

}; //namespace itensor
