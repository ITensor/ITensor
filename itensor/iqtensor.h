#ifndef __IQ_H
#define __IQ_H
#include "itensor.h"
#include "iqindex.h"
#include <set>
#include <list>
#include <map>

class IQTDat
{
    typedef std::list<ITensor>::iterator       iten_it;
    typedef std::list<ITensor>::const_iterator const_iten_it;
private:
    mutable unsigned int numref;
    mutable bool rmap_init;
public:
    mutable std::list<ITensor> itensor; // This is mutable to allow reordering
    std::vector<IQIndex> iqindex_;
    mutable std::map<ApproxReal,iten_it> rmap; //mutable so that const IQTensor methods can use rmap

    IQTDat() : numref(0), rmap_init(false) { }

    explicit IQTDat(const IQIndex& i1)
    : numref(0), rmap_init(false), iqindex_(1)
    { 
        iqindex_[0] = i1; 
    }

    IQTDat(const IQIndex& i1, const IQIndex& i2)
    : numref(0), rmap_init(false), iqindex_(2)
    { 
        iqindex_[0] = i1; 
        iqindex_[1] = i2; 
    }

    IQTDat(const IQIndex& i1, const IQIndex& i2, const IQIndex& i3)
    : numref(0), rmap_init(false), iqindex_(3)
    { 
        iqindex_[0] = i1; 
        iqindex_[1] = i2; 
        iqindex_[2] = i3; 
    }

    IQTDat(const IQIndex& i1, const IQIndex& i2, 
           const IQIndex& i3, const IQIndex& i4,
           const IQIndex& i5 = IQIndex::Null(), 
           const IQIndex& i6 = IQIndex::Null(), 
           const IQIndex& i7 = IQIndex::Null(), 
           const IQIndex& i8 = IQIndex::Null())
    : numref(0), rmap_init(false), iqindex_(4)
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

    explicit IQTDat(std::vector<IQIndex>& iqinds_) : numref(0), rmap_init(false) { iqindex_.swap(iqinds_); }

    explicit IQTDat(const IQTDat& other) 
    : numref(0), rmap_init(false), itensor(other.itensor), iqindex_(other.iqindex_)
    { }

    explicit IQTDat(std::istream& s) : numref(0) { read(s); }

    void read(std::istream& s)
    {
        uninit_rmap();
        size_t size;
        s.read((char*) &size,sizeof(size));
        itensor.resize(size);
        foreach(ITensor& t, itensor) t.read(s);

        s.read((char*) &size,sizeof(size));
        iqindex_.resize(size);
        foreach(IQIndex& I, iqindex_) I.read(s);
    }

    void write(std::ostream& s) const
    {
        size_t size = itensor.size();
        s.write((char*) &size,sizeof(size));
        foreach(const ITensor& t, itensor) t.write(s);

        size = iqindex_.size();
        s.write((char*) &size,sizeof(size));
        foreach(const IQIndex& I, iqindex_) I.write(s);
    }

    void init_rmap() const
    {
        if(rmap_init) return;
        for(iten_it it = itensor.begin(); it != itensor.end(); ++it)
        { rmap[ApproxReal(it->unique_Real())] = it; }
        rmap_init = true;
    }

    void uninit_rmap() const 
    { 
        assert(numref <= 1); 
        rmap.clear();
        rmap_init = false; 
    }

    bool has_itensor(const ApproxReal& r) const
    { 
        init_rmap();
        return rmap.count(r) == 1; 
    }

    void insert_itensor(const ApproxReal& r, const ITensor& t)
    {
        itensor.push_front(t);
        rmap[r] = itensor.begin();
    }

    void clean(Real min_norm)
    {
        std::list<ITensor> nitensor;
        foreach(const ITensor& t, itensor)
        {
            if(t.norm() >= min_norm)
                nitensor.push_back(t);
            /*
            else
            {
                std::cout << boost::format("Discarding tensor with norm %.2E")%t.norm() << std::endl;
                PrintDat(t);
            }
            */
        }
        itensor.swap(nitensor);
    }

    inline void* operator new(size_t size) throw(std::bad_alloc)
        { return allocator.alloc(); }

    inline void operator delete(void* p) throw()
        { return allocator.dealloc(p); }

    ENABLE_INTRUSIVE_PTR(IQTDat)
private:
    static DatAllocator<IQTDat> allocator;
    ~IQTDat() { } //must be dynamically allocated
    void operator=(const IQTDat&);
};

class IQCombiner;

class IQTensor
{
public:
    typedef IQIndex IndexT;
    typedef IQIndexVal IndexValT;
    typedef IQCombiner CombinerT;
    typedef std::list<ITensor>::iterator iten_it;
    typedef std::list<ITensor>::const_iterator const_iten_it;
    typedef std::vector<IQIndex>::iterator iqind_it;
    typedef std::vector<IQIndex>::const_iterator const_iqind_it;
    static const IQIndex& ReImIndex()
    {
        return IQIndex::IndReIm();
    }
private:
    boost::intrusive_ptr<IQTDat> p;

    void solo()
    {
        assert(p != 0);
        if(p->count() != 1)
        {
            boost::intrusive_ptr<IQTDat> new_p(new IQTDat(*p));
            p.swap(new_p);
        }
    }
public:
    int r() const { return p->iqindex_.size(); }
    inline const IQIndex& index(int j) const { return GET(p->iqindex_,j-1); }
    inline int iten_size() const { return p->itensor.size(); }
    inline bool iten_empty() const { return p->itensor.empty(); }
    inline bool is_null() const { return p == 0; }
    inline bool is_not_null() const { return p != 0; }
    inline int num_index() const { return p->iqindex_.size(); }
    
    //----------------------------------------------------
    //IQTensor: iterators 
    const_iten_it const_iten_begin() const { return p->itensor.begin(); }
    const_iten_it const_iten_end() const { return p->itensor.end(); }
    std::pair<const_iten_it,const_iten_it> itensors() const { return std::make_pair(p->itensor.begin(),p->itensor.end()); }

    const_iqind_it const_iqind_begin() const { return p->iqindex_.begin(); }
    const_iqind_it const_iqind_end()   const { return p->iqindex_.end(); }
    std::pair<const_iqind_it,const_iqind_it> iqinds() const { return std::make_pair(p->iqindex_.begin(),p->iqindex_.end()); }

    //----------------------------------------------------
    //IQTensor: Constructors

    IQTensor() : p(0) {}

    explicit IQTensor(const IQIndex& i1) 
    : p(new IQTDat(i1))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2) 
    : p(new IQTDat(i1,i2))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3) 
    : p(new IQTDat(i1,i2,i3))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4) 
    : p(new IQTDat(i1,i2,i3,i4))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5)
    : p(new IQTDat(i1,i2,i3,i4,i5))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5,const IQIndex& i6)
    : p(new IQTDat(i1,i2,i3,i4,i5,i6))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
             const IQIndex& i7)
    : p(new IQTDat(i1,i2,i3,i4,i5,i6,i7))
    { }

    IQTensor(const IQIndex& i1,const IQIndex& i2,const IQIndex& i3,
             const IQIndex& i4,const IQIndex& i5,const IQIndex& i6,
             const IQIndex& i7,const IQIndex& i8)
    : p(new IQTDat(i1,i2,i3,i4,i5,i6,i7,i8))
    { }

    explicit IQTensor(std::vector<IQIndex>& iqinds_) 
    : p(new IQTDat(iqinds_))
    { }

    IQTensor(const IQIndexVal& iv1) 
    : p(new IQTDat(iv1.iqind))
    { 
        operator()(iv1) = 1;
    }

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2) 
    : p(new IQTDat(iv1.iqind,iv2.iqind))
    { 
        operator()(iv1,iv2) = 1;
    }

    IQTensor(const IQIndexVal& iv1, const IQIndexVal& iv2,
             const IQIndexVal& iv3) 
    : p(new IQTDat(iv1.iqind,iv2.iqind,iv3.iqind))
    { 
        operator()(iv1,iv2,iv3) = 1;
    }

    IQTensor(ITmaker itm) : p(new IQTDat(IQIndex::IndReIm()))
    {
        if(itm == makeComplex_1)      
            operator+=(ITensor::Complex_1());
        else if(itm == makeComplex_i) 
            operator+=(ITensor::Complex_i());
    }

    IQTensor(IQmaker i) : p(new IQTDat())
    {
        Index s("sing");
        IQIndex single("single",s,QN());
        p->iqindex_.push_back(single);
        ITensor st(s,1);
        operator+=(st);
    }

    IQTensor(PrimeType pt,const IQTensor& other) : p(other.p)
    { doprime(pt); }

    explicit IQTensor(std::istream& s) : p(0) { read(s); }

    static const IQTensor& Sing()
    {
        static const IQTensor Sing_(makeSing);
        return Sing_;
    }

    static const IQTensor& Complex_1()
    {
        static const IQTensor Complex_1_(makeComplex_1);
        return Complex_1_;
    }

    static const IQTensor& Complex_i()
    {
        static const IQTensor Complex_i_(makeComplex_i);
        return Complex_i_;
    }

    void read(std::istream& s)
    {
        bool null_;
        s.read((char*) &null_,sizeof(null_));
        if(null_) { *this = IQTensor(); return; }
        p = new IQTDat(s);
    }

    void write(std::ostream& s) const
    {
        bool null_ = is_null();
        s.write((char*) &null_,sizeof(null_));
        if(null_) return;
        p->write(s);
    }


    //----------------------------------------------------
    //IQTensor operators

    IQTensor operator*(IQTensor other) const { other *= *this; return other; }
    IQTensor& operator*=(const IQTensor& other);

    IQTensor& operator+=(const IQTensor& o);
    IQTensor operator+(const IQTensor& o) const { IQTensor res(*this); res += o; return res; }
    IQTensor operator-(const IQTensor& o) const 
    { IQTensor res(o); res *= -1; res += *this; return res; }

    IQTensor& operator*=(Real fac) 
    { 
        solo();
        if(fac == 0) { p->itensor.clear(); p->uninit_rmap(); return *this; }
        foreach(ITensor& t, p->itensor) { t *= fac; }
        return *this; 
    }
    IQTensor operator*(Real fac) { IQTensor res(*this); res *= fac; return res; }
    friend inline IQTensor operator*(Real fac, IQTensor T) { T *= fac; return T; }

    /*
    operator ITensor() const
    {
        std::vector<Index> indices;
        foreach(const IQIndex& I, p->iqindex_)
        {
            if(I.type() != Site) 
            { Error("IQTensor to ITensor conversion requires all IQIndex's of type Site."); }
            indices.push_back(I);
        }
        ITensor res(indices);
        match_order();
        return res;
    }
    */

    void insert(const ITensor& t) 
    { 
        solo();
        ApproxReal r(t.unique_Real());
        if(p->has_itensor(r))
        {
            Print(*(p->rmap[r])); Print(t);
            Error("Can't insert ITensor with identical structure twice, use operator+=.");
        }
        p->insert_itensor(r,t);
    } //IQTensor::insert(const ITensor& t)

    IQTensor& operator+=(const ITensor& t) 
    { 
        solo();
        ApproxReal r(t.unique_Real());

        if(t.scale().isRealZero()) { return *this; }

        if(!p->has_itensor(r)) { p->insert_itensor(r,t); }
        else { *(p->rmap[r]) += t; }
        return *this;
    } //IQTensor::operator+=(ITensor)

    Real& operator()(const IQIndexVal& iv1, const IQIndexVal& iv2 = IQIndexVal::Null(), const IQIndexVal& iv3 = IQIndexVal::Null(),
                     const IQIndexVal& iv4 = IQIndexVal::Null(), const IQIndexVal& iv5 = IQIndexVal::Null(), const IQIndexVal& iv6 = IQIndexVal::Null(),
                     const IQIndexVal& iv7 = IQIndexVal::Null(), const IQIndexVal& iv8 = IQIndexVal::Null())
	{
        solo();
        boost::array<IQIndexVal,NMAX+1> iv 
            = {{ IQIndexVal::Null(), iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        Real ur = 0; int nn = 0; 
        while(GET(iv,nn+1).iqind != IQIndexVal::Null().iqind) 
        { ++nn; ur += GET(iv,nn).index().unique_Real(); }
        if(nn != r()) Error("Not enough IQIndexVals provided");
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
        return (p->rmap[r])->operator()(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8);
    }

    //----------------------------------------------------
    //IQTensor quantum number methods

    QN div() const
    {
        QN div_;
        assert(p != 0);
        if(p->itensor.empty())
        {   
            this->printIQInds("this");
            Error("IQTensor has no blocks");
        }
        const ITensor& t = p->itensor.front();
        for(int j = 1; j <= t.r(); ++j)
            { div_ += qn(t.index(j))*dir(t.index(j)); }
        return div_;
    }

    bool checkDiv(QN expected) const
    {
        assert(p != 0);
        if(p->itensor.empty())
        {   
            this->printIQInds("this");
            Error("IQTensor has no blocks");
        }
        foreach(const ITensor& t, p->itensor)
        {
            QN div_;
            for(int j = 1; j <= t.r(); ++j)
                { div_ += qn(t.index(j))*dir(t.index(j)); }
            if(div_ != expected)
            {
                std::cerr << "Block didn't match expected div\n";
                Print(t);
                return false;
            }
        }
        return true;
    }

    QN qn(const Index& in) const
    {
        int iqq = find_iqind(in)-1;
        if(iqq == -1) Error("qn: cant find index");
        return p->iqindex_[iqq].qn(in);
    } //end IQTensor::qn

    //void negate_qn() { foreach(IQIndex& I, p->iqindex_) I.negate(); }

    Arrow dir(const Index& in) const
    {
        int iqq = find_iqind(in)-1;
        if(iqq == -1) 
        {
            this->print("this IQTensor");
            in.print("in"); 
            Error("IQTensor::dir(Index&): cant find Index in IQIndices");
        }
        return p->iqindex_[iqq].dir();
    } //end IQTensor::dir


    //----------------------------------------------------
    //IQTensor: prime methods

    void ind_inc_prime(const IQIndex& i,int inc)
    {
        solo();
        p->uninit_rmap();
        bool gotit = false;
        foreach(IQIndex& jj, p->iqindex_)
        if(jj.noprime_equals(i))
        {
            gotit = true;
            int p = jj.primeLevel();
            jj.mapprime(p,p+inc);
        }

        if(!gotit)
        {
            std::cerr << "IQIndex was " << i << "\n";
            Error("ind_inc_prime: couldn't find IQIndex");
        }

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
        for(int ii = 1; ii <= jj->r(); ++ii)
        if(i.hasindex_noprime(jj->index(ii)))
        {
            int p = jj->index(ii).primeLevel();
            jj->mapprimeind(jj->index(ii),p,p+inc);
        }
    } //end IQTensor::ind_inc_prime

    void noprime(PrimeType pt = primeBoth)
    {
        solo();
        p->uninit_rmap();
        foreach(IQIndex& J, p->iqindex_)
	    J.noprime(pt);

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    foreach(ITensor& t, p->itensor) t.noprime(pt);
    } //end IQTensor::noprime
    friend inline IQTensor deprimed(IQTensor A)
    { A.noprime(); return A; }

    void noprimelink()
    {
        solo();
        p->uninit_rmap();
        foreach(IQIndex& J, p->iqindex_)
	    if(J.type() == Link) J.noprime();

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
	    foreach(ITensor& t, p->itensor) t.noprime(primeLink);
    } //end IQTensor::noprimelink

    void doprime(PrimeType pt)
    {
        solo();
        p->uninit_rmap();

        DoPrimer prim(pt);
        for_each(p->iqindex_.begin(), p->iqindex_.end(),prim);

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
            { jj->doprime(pt); }
    } //end IQTensor::doprime

    //no need to keep prime level small
    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
    {
        solo();
        p->uninit_rmap();

        MapPrimer prim(plevold,plevnew,pt);
        for_each(p->iqindex_.begin(), p->iqindex_.end(),prim);

        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
            { jj->mapprime(plevold,plevnew,pt); }
    } //end IQTensor::mapprime

    void primeind(const IQIndex& I)
    {
        solo();
        p->uninit_rmap();
        foreach(IQIndex& J, p->iqindex_)
        { if(J == I) J = J.primed(); }

        foreach(ITensor& t, p->itensor)
        foreach(const inqn& x, I.iq())
        { if(t.hasindex(x.index)) t.primeind(x.index); }
    }
    friend inline IQTensor primeind(IQTensor A, const IQIndex& I)
    { A.primeind(I); return A; }
    friend inline IQTensor primeind(IQTensor A, const IQIndex& I, 
                                                const IQIndex& J)
    { A.primeind(I); A.primeind(J); return A; }

    friend inline IQTensor primed(IQTensor A)
    { A.doprime(primeBoth); return A; }

    void primesite() { doprime(primeSite); }
    friend inline IQTensor primesite(IQTensor A)
    { A.doprime(primeSite); return A; }

    void primelink() { doprime(primeLink); }
    friend inline IQTensor primelink(const IQTensor& A)
    { IQTensor res(A); res.doprime(primeLink); return res; }


    //----------------------------------------------------
    //IQTensor index methods

    int find_iqind(const Index& I) const
    {
        for(size_t j = 0; j < p->iqindex_.size(); ++j)
        if(p->iqindex_[j].hasindex(I)) { return j+1; }
        return 0;
    }

    //Return true if one of the ITensors uses this Index
    bool uses_ind(const Index& i) const
    {
        foreach(const ITensor& it, p->itensor)
        if(it.hasindex(i)) { return true; }
        return false;
    }

    int findindex(const IQIndex& i) const
    {
        std::vector<IQIndex>::const_iterator f = find(p->iqindex_.begin(),p->iqindex_.end(),i);
        if(f == p->iqindex_.end()) return 0;
        else return (f - p->iqindex_.begin())+1;
    }

    bool hastype(IndexType t) const
    {
        foreach(const IQIndex& I, p->iqindex_)
        if(I.type() == t) { return true; }
        return false;
    }

    const IQIndex& findtype(IndexType t) const
    {
        foreach(const IQIndex& I, p->iqindex_)
        if(I.type() == t) { return I; }
        Error("IQTensor::findtype: couldn't find type");
        return IQIndex::Null();
    }

    const IQIndex& finddir(Arrow dir) const
    {
        foreach(const IQIndex& I, p->iqindex_)
        if(I.dir() == dir) { return I; }
        Error("IQTensor::finddir: couldn't find dir");
        return IQIndex::Null();
    }

    bool hasindex(const IQIndex& i) const 
    { 
        for(size_t j = 0; j < p->iqindex_.size(); ++j)
            { if(i == p->iqindex_[j]) return true; }
        return false;
    }

    bool is_complex() const { return findindex(IQIndex::IndReIm()) != 0; }

    void addindex1(const IQIndex& I)
    {
        if(I.m() != 1) Error("IQTensor::operator*=(IQIndex): IQIndex must have m == 1.");    
        solo(); p->uninit_rmap();
        foreach(ITensor& t, p->itensor) t.addindex1(I.index(1));
        p->iqindex_.push_back(I);
    }

    //----------------------------------------------------
    //IQTensor miscellaneous methods

    Real unique_Real() const
    {
        Real ur = 0.0;
        foreach(const IQIndex& I, p->iqindex_) ur += I.unique_Real();
        return ur;
    }
    int num_index(IndexType t) const 
    { 
        int count = 0;
        foreach(const IQIndex& I, p->iqindex_)
        if(I.type() == t) ++count;
        return count;
    }

    Real norm() const
    {
        Real res = 0;
        foreach(const ITensor& t, p->itensor) 
            res += sqr(t.norm());
        return sqrt(res);
    }

    Real sumels() const
    {
        Real res = 0;
        foreach(const ITensor& t, p->itensor)
            res += t.sumels();
        return res;
    }

    void scaleOutNorm() const
    {
        Real f = norm();
        foreach(const ITensor& t, p->itensor)
            { t.scaleTo(f); }
    }

    void scaleTo(LogNumber newscale) const
    {
        foreach(const ITensor& t, p->itensor)
            { t.scaleTo(newscale); }
    }

    inline void clean(Real min_norm = MIN_CUT)
    { 
        solo();
        p->clean(min_norm); 
    }

    int vec_size() const
    {
        int s = 0;
        for(const_iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
            s += jj->vec_size();
        return s;
    }
    void assignToVec(VectorRef v) const
    {
        if(vec_size() != v.Length())
            Error("Mismatched sizes in IQTensor::assignToVec(VectorRef v).");
        int off = 1;
        for(const_iten_it jj = const_iten_begin(); jj != const_iten_end(); ++jj)
        {
            int d = jj->vec_size();
            jj->assignToVec(v.SubVector(off,off+d-1));
            off += d;
        }
    }
    void assignFromVec(VectorRef v)
    {
        solo();
        if(vec_size() != v.Length())
            Error("bad size");
        int off = 1;
        for(iten_it jj = p->itensor.begin(); jj != p->itensor.end(); ++jj)
        {
            int d = jj->vec_size();
            jj->assignFromVec(v.SubVector(off,off+d-1));
            off += d;
        }
    }
    void GetSingComplex(Real& re, Real& im) const;

    
    void Randomize() { solo(); foreach(ITensor& t, p->itensor) t.Randomize(); }

    void print(std::string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); std::cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    void printIQInds(std::string name = "") const
    { 
        std::cerr << "\n" << name << " (IQIndices only) = \n";
        for(size_t j = 0; j < p->iqindex_.size(); ++j)
        { std::cerr << p->iqindex_[j] << "\n\n"; }
        std::cerr << "---------------------------\n\n";
    }

    void assignFrom(const IQTensor& other) const
    {
        std::map<ApproxReal,iten_it> semap;
        for(iten_it i = p->itensor.begin(); i != p->itensor.end(); ++i)
        { semap[ApproxReal(i->unique_Real())] = i; }
        for(const_iten_it i = other.p->itensor.begin(); i != other.p->itensor.end(); ++i)
        {
            ApproxReal se = ApproxReal(i->unique_Real());
            if(semap.count(se) == 0)
            {
                std::cout << "warning assignFrom semap.count is 0" << std::endl;
                std::cerr << "offending ITensor is " << *i << "\n";
                Error("bad assignFrom count se");
            }
            else { semap[se]->assignFrom(*i); }
        }
    }

    void SplitReIm(IQTensor& re, IQTensor& im) const;

    void conj()
    {
        solo();
        if(!is_complex())
        { foreach(IQIndex& I,p->iqindex_) I.conj(); }
        else
        {
            IQTensor r,i;
            SplitReIm(r,i);
            foreach(IQIndex& I,r.p->iqindex_) I.conj();
            foreach(IQIndex& I,i.p->iqindex_) I.conj();
            i *= -1.0;
            *this = r * IQTensor::Complex_1() + IQTensor::Complex_i() * i;
        }
    }

private:
    void match_order() const
    {
        const int s = p->iqindex_.size();
        foreach(ITensor& t, p->itensor)
        {
            if(t.r() != s)
            {
                Print(*this); Print(t);
                Error("match_order: ds not same");
            }
            Permutation P;
            for(int j = 1; j <= t.r(); ++j)
            {
                bool gotone = false;
                for(int k = 0; k < s; ++k)
                if(p->iqindex_[k].hasindex(t.index(j)))
                { P.from_to(j,k+1); gotone = true; break; }
                if(!gotone) Error("match_order: !gotone");
            }
            t.reshape(P);
        }
    }

public:

    inline friend std::ostream& operator<<(std::ostream & s, const IQTensor &t)
    {
        s << "\n----- IQTensor -----\nIQIndices: " << std::endl;
        foreach(const IQIndex& I, t.iqinds()) s << "  " << I << std::endl;
        s << "ITensors: \n";
        foreach(const ITensor& i, t.itensors())
        s <<"	" << i << "\n";
        s << "-------------------" << "\n\n";
        return s;
    }

}; //class IQTensor

inline Real ReSingVal(const IQTensor& x)
{
    Real re, im;
    x.GetSingComplex(re,im);
    return re;
}
inline Real Dot(const IQTensor& x, const IQTensor& y, bool doconj = true)
{
    IQTensor res(IQTensor::Sing()*(doconj ? conj(x) : x)*y);
    return ReSingVal(res);
}
inline void Dot(const IQTensor& x, const IQTensor& y, Real& re, Real& im, bool doconj = true)
{
    IQTensor res(IQTensor::Sing()*(doconj ? conj(x) : x)*y);
    res.GetSingComplex(re,im);
}

inline bool checkQNs(const ITensor& t) { return true; }
//Checks if all IQTensor blocks have the same divergence
inline bool checkQNs(const IQTensor& T)
{
    QN qtot = T.div();
    foreach(const ITensor& it, T.itensors())
    {
        QN q;
        for(int j = 1; j <= it.r(); ++j) 
            { q += T.qn(it.index(j))*T.dir(it.index(j)); }

        if(q != qtot) 
        {
            std::cerr << "checkQNs: inconsistent QN.\n";
            std::cerr << "\nqtot = " << qtot << "\n\n";
            std::cerr << "Offending ITensor = " << it << "\n\n";
            T.printIQInds("T");
            return false;
        }
    }
    return true;
}

template<class T> class Printit
{
public:
    std::ostream& s;
    std::string spacer;
    Printit(std::ostream& _s, std::string _spacer) : s(_s), spacer(_spacer) {}
    void operator()(const T& t) { s << t << spacer; }
};


class SiteOp
{
    const IQIndex si;
    mutable bool made_iqt, made_t;
    mutable IQTensor iqt;
    mutable ITensor t;
    std::map<ApproxReal, std::pair<IQIndexVal,IQIndexVal> > ivmap;
    std::map<ApproxReal, Real> valmap;
    typedef std::map<ApproxReal, Real>::value_type valmap_vt;

    std::pair<IQIndexVal,IQIndexVal> civmap(const ApproxReal& key) const { return ivmap.find(key)->second; }

    void make_iqt() const
    {
        if(made_iqt) return;
        iqt = IQTensor(conj(si),si.primed());
        foreach(const valmap_vt& x, valmap)
        { iqt(civmap(x.first).first,civmap(x.first).second) = x.second; }
        made_iqt = true;
    }
    void make_t() const
    {
        if(made_t) return;
        t = ITensor(conj(si),si.primed());
        foreach(const valmap_vt& x, valmap)
        { t(civmap(x.first).first,civmap(x.first).second) = x.second; }
        made_t = true;
    }
public:
    SiteOp(const IQIndex& si_) : si(si_), made_iqt(false), made_t(false) { }
    operator ITensor() const { make_t(); return t; }
    operator IQTensor() const { make_iqt(); return iqt; }

    Real& operator()(const IQIndexVal& iv, const IQIndexVal& ivp)
    {
        if(iv.iqind.primeLevel() != 0) Error("SiteOp::operator(): first IndexVal must be unprimed.");
        if(ivp.iqind.primeLevel() != 1) Error("SiteOp::operator(): second IndexVal must be primed.");
        Real r = iv.iqind.unique_Real() + ivp.iqind.unique_Real() + iv.i + 1000*ivp.i;
        ivmap[ApproxReal(r)] = std::make_pair(iv,ivp);
        return valmap[ApproxReal(r)];
    }

    void print(std::string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); std::cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }
    friend std::ostream& operator<<(std::ostream& s, const SiteOp& op) { return s << IQTensor(op); }

    template <class X>
    ITensor operator*(const X& x) const { ITensor res(*this); res *= x; return res; }
    template <class X>
    friend inline ITensor operator*(const X& x, const SiteOp& op) { ITensor res(op); return (res *= x); }

    IQTensor operator*(const IQTensor& t) const { IQTensor res(*this); res *= t; return res; }
    friend inline IQTensor operator*(const IQTensor& t, const SiteOp& op) { IQTensor res(op); return (res *= t); }

};

#endif
