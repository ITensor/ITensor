#include "iqtsparse.h"
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


//
// IQTSparse
//

IQTSparse::
IQTSparse()
    :
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
	assert(d_ != 0);
	if(d_->count() != 1)
	    {
	    boost::intrusive_ptr<IQTSDat> new_d_(new IQTSDat(*d_));
	    d_.swap(new_d_);
	    }
    }

void IQTSparse::
soloIndex()
    {
	assert(is_ != 0);
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
    if(!S.isDiag()) 
        Error("product only implemented for diagonal IQTSparse's");

    } // void product(const IQTSparse& S, const ITensor& T, ITensor& res)

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
