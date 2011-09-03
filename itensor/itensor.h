#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "types.h"
#include "real.h"
#include "index.h"
#include "permutation.h"
#include ".profiling/prodstats.h"
#include ".profiling/count_copies.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using std::stringstream;
using std::pair;
using std::make_pair;
using std::vector;
using boost::format;
using boost::array;

enum ITmaker {makeComplex_1,makeComplex_i,makeConjTensor};

class Permutation;

class Counter
{
private:
    void reset(int a)
	{
        i.assign(a);
        ind = 1;
	}
    int rn_,r_;
public:
    array<int,NMAX+1> n;
    array<int,NMAX+1> i;
    int ind;

    Counter() : rn_(0)
	{
        n.assign(1); n[0] = 0;
        reset(0);
	}

    Counter(const array<Index,NMAX+1>& ii,int rn,int r) { init(ii,rn,r); }

    void init(const array<Index,NMAX+1>& ii, int rn, int r)
    {
        rn_ = rn;
        r_ = r;
        n[0] = 0;
        for(int j = 1; j <= rn_; ++j) 
        { GET(n,j) = ii[j].m(); }
        for(int j = rn_+1; j <= NMAX; ++j) 
        { n[j] = 1; }
        reset(1);
    }

    Counter& operator++()
	{
        ++ind;
        ++i[1];
        if(i[1] > n[1])
        for(int j = 2; j <= rn_; ++j)
        {
            i[j-1] = 1;
            ++i[j];
            if(i[j] <= n[j]) break;
        }
        if(i[rn_] > n[rn_]) reset(0); //set 'done' condition
        return *this;
	}

    bool operator!=(const Counter& other) const
	{
        for(int j = 1; j <= NMAX; ++j)
        { if(i[j] != other.i[j]) return true; }
        return false;
	}
    bool operator==(const Counter& other) const
	{ return !(*this != other); }

    //static const Counter *pend;
    //static const Counter& done;
    static const Counter done;

    friend inline ostream& operator<<(ostream& s, const Counter& c)
    {
        s << "("; for(int i = 1; i < c.r_; ++i){s << c.i[i] << " ";} s << c.i[c.r_] << ")";
        return s;
    }
};
#ifdef THIS_IS_MAIN
//const Counter* Counter::pend = new Counter;
//const Counter& Counter::done = *pend;
const Counter Counter::done;
#endif                                  

//#define DO_ALT
#ifdef DO_ALT
struct PDat
{
    Permutation I; 
    Vector v;
    PDat(const Permutation& P_, const Vector& v_) : I(P_.inverse()), v(v_) { }
    PDat(const Permutation& P_) : I(P_.inverse()) { }
};
#endif

//Storage for ITensors
class ITDat
{
private:
    mutable unsigned int numref;
public:
    Vector v;
#ifdef DO_ALT
    vector<PDat> alt;
#endif

    ITDat() : numref(0), v(0) { }

    explicit ITDat(int size) 
    : numref(0), v(size)
	{ assert(size > 0); v = 0; }

    explicit ITDat(const Vector& v_) 
    : numref(0), v(v_)
    { }

    explicit ITDat(Real r) 
    : numref(0), v(1)
    { v = r; }

    explicit ITDat(istream& s) : numref(0) { read(s); }

    explicit ITDat(const ITDat& other) 
    : numref(0), v(other.v)
    { }

    void read(istream& s)
	{ 
        int size = 0;
        s.read((char*) &size,sizeof(size));
        v.ReDimension(size);
        s.read((char*) v.Store(), sizeof(Real)*size);
    }

    void write(ostream& s) const 
    { 
        const int size = v.Length();
        s.write((char*) &size, sizeof(size));
        s.write((char*) v.Store(), sizeof(Real)*size); 
    }

    void print() const { cout << "ITDat: v = " << v; }

    friend class ITensor;
    ENABLE_INTRUSIVE_PTR(ITDat)
private:
    void operator=(const ITDat&);
    ~ITDat() { } //must be dynamically allocated
};

class Combiner;

class ITensor; extern ITensor Complex_1, Complex_i, ConjTensor;

class ITensor
{
public:
    typedef Index IndexT;
    typedef IndexVal IndexValT;
    typedef Combiner CombinerT;
    typedef array<Index,NMAX+1>::const_iterator indexn_it;
    static const Index& ReImIndex;
private:
    mutable intrusive_ptr<ITDat> p; //mutable: const methods may want to reshape data
    int rn;
    mutable array<Index,NMAX+1> _indexn; //Indices having m!=1, maximum of 8 (_indexn[0] not used), mutable to allow reordering
    mutable vector<Index>       _index1; //Indices having m==1
    Real ur;
    mutable LogNumber scale_; //mutable since e.g. normlogto is logically const

    void allocate(int dim) { p = new ITDat(dim); }
    void allocate() { p = new ITDat(); }

#ifdef DO_ALT
    void newAltDat(const Permutation& P) const
    { p->alt.push_back(PDat(P)); }
    PDat& lastAlt() const { return p->alt.back(); } 
#endif

    //Disattach self from current ITDat and create own copy instead.
    //Necessary because ITensors logically represent distinct
    //objects even though they may share data in reality.
    void solo() const
	{
        assert(p != 0);
        if(p->count() != 1) 
        {
            intrusive_ptr<ITDat> new_p(new ITDat(*p));
            p.swap(new_p);
            IF_COUNT_COPIES(++copycount;)
        }
	}
    
    void set_unique_Real()
	{
        ur = 0.0;
        for(int j = 1; j <= rn; ++j)
        { ur += GET(_indexn,j).unique_Real(); }
        foreach(const Index& I, _index1)
        { ur += I.unique_Real(); }
	}

    void _construct2(const Index& i1, const Index& i2)
    {
        if(i1.m()==1) _index1.push_back(i1); else { GET(_indexn,++rn) = i1; }
        if(i2.m()==1) _index1.push_back(i2); else { GET(_indexn,++rn) = i2; }
        allocate(i1.m()*i2.m()); 
        set_unique_Real();
    }

    //Prefer to map via a Combiner
    //'mapindex' is useful and efficient if used in a safe way, however.
    void mapindex(const Index& i1, const Index& i2)
	{
        assert(i1.m() == i2.m());
        for(int j = 1; j <= rn; ++j) 
        if(GET(_indexn,j) == i1) 
        {
            solo();
            GET(_indexn,j) = i2;
            set_unique_Real();
            return;
        }
        vector<Index>::iterator i1pos = find(_index1.begin(),_index1.end(),i1);
        if(i1pos == _index1.end()) 
        {
            cerr << "\nFor ITensor = " << *this << "\n";
            cerr << "Missing index i1 = " << i1 << "\n";
            Error("ITensor::mapindex(i1,i2): ITensor does not have index i1");
        }
        *i1pos = i2;
        set_unique_Real();
	}

    void getperm(const array<Index,NMAX+1>& oth_indexn, Permutation& P) const
    {
        for(int j = 1; j <= rn; ++j)
        {
            bool got_one = false;
            for(int k = 1; k <= rn; ++k)
            if(oth_indexn[j] == _indexn[k])
            { P.from_to(j,k); got_one = true; break; }
            if(!got_one)
            {
                cerr << "j = " << j << "\n";
                Print(*this);
                foreach(const Index& I, oth_indexn) { cerr << I << "\n"; }
                Error("getpermBtoA: no matching index");
            }
        }
    }

    friend void toMatrixProd(const ITensor& L, const ITensor& R, 
                             array<bool,NMAX+1>& contractedL, array<bool,NMAX+1>& contractedR, 
                             MatrixRefNoLink& lref, MatrixRefNoLink& rref);

    int fillFromIndices(const vector<Index>& I)
    {
        assert(rn == 0);
        assert(_index1.empty());
        int alloc_size = 1;
        for(size_t n = 0; n < I.size(); ++n)
        {
            const Index& i = I[n];
            if(i == IndNull) Error("ITensor: null Index in constructor.");
            if(i.m()==1) _index1.push_back(i); 
            else 
            { 
                DO_IF_DEBUG(if(rn == NMAX) Error("ITensor(const vector<Index>& I): too many indices with m > 1");)
                GET(_indexn,++rn) = i; 
                alloc_size *= i.m(); 
            }
        }
        return alloc_size;
    }

public:

    //Accessor Methods ----------------------------------------------

    Real unique_Real() const { return ur; }	// depends on indices only, unordered
    const Index& index(int j) const { return (j > rn ? GET(_index1,j-(rn+1)) : GET(_indexn,j)); }
    const Index& indexn(int j) const { return GET(_indexn,j); }
    int r() const { return rn + _index1.size(); }
    int r_n() const { return rn; }
    int r_1() const { return _index1.size(); }
    inline int m(int j) const { return (j > rn ? 1 : GET(_indexn,j).m()); }

    bool is_null() const { return (p == 0); }
    bool is_not_null() const { return (p != 0); }
    bool is_complex() const { return findindexn(IndReIm) > 0; }
    bool is_not_complex() const { return (findindexn(IndReIm) == 0); }

    LogNumber scale() const { return scale_; }

    //These methods can be used for const iteration over Indices in a foreach loop
    //e.g. foreach(const Index& I, t.indexn() ) { ... }
    const pair<indexn_it,indexn_it> indexn() const { return make_pair(_indexn.begin()+1,_indexn.begin()+rn+1); }
    const vector<Index>&            index1() const { return _index1; }

    void initCounter(Counter& C) const { C.init(_indexn,rn,r()); }

    //Constructors --------------------------------------------------

    ITensor() : p(0), rn(0), ur(0)  { }

    ITensor(Real val) : rn(0)
	{ 
        allocate(1);
        p->v = val;
        set_unique_Real();
    }

    explicit ITensor(const Index& i1) : rn(0)
	{ 
        if(i1.m()==1) _index1.push_back(i1); else { _indexn[1] = i1; ++rn; }
        allocate(i1.m());
        set_unique_Real();
    }

    ITensor(const Index& i1, Real val) : rn(0)
	{ 
        if(i1.m()==1) _index1.push_back(i1); else { _indexn[1] = i1; ++rn; }
        allocate(i1.m()); p->v = val; 
        set_unique_Real();
    }

    ITensor(const Index& i1, const Vector& V) : p(new ITDat(V)), rn(0)
	{ 
        if(i1.m() != V.Length()) Error("Mismatch of Index and Vector sizes.");
        if(i1.m()==1) _index1.push_back(i1); else { GET(_indexn,1) = i1; ++rn; }
        set_unique_Real();
    }

    ITensor(Index i1,Index i2) : rn(0)
	{ _construct2(i1,i2); }

    //Create an ITensor as a matrix with 'a' on the diagonal
    ITensor(Index i1,Index i2,Real a) : rn(0)
    {
        _construct2(i1,i2);
        int nn = min(i1.m(),i2.m());
        if(rn == 2) for(int i = 1; i <= nn; ++i) p->v((i-1)*i1.m()+i) = a;
        else p->v(1) = a;
    }

    ITensor(Index i1,Index i2,const MatrixRef& M) : rn(0)
    {
        _construct2(i1,i2);
        if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
        { Error("ITensor(Index,Index,Matrix): Mismatch of Index sizes and matrix."); }
        MatrixRef dref; p->v.TreatAsMatrix(dref,i2.m(),i1.m());
        dref = M.t();
    }

    ITensor(Index i1, Index i2, Index i3,
            Index i4 = IndNull, Index i5 = IndNull, Index i6 = IndNull,
            Index i7 = IndNull, Index i8 = IndNull)
            : rn(0)
    {
        array<Index,NMAX+1> ii = {{ IndNull, i1, i2, i3, i4, i5, i6, i7, i8 }};
        int dim = 1;
        for(int n = 1; n <= NMAX; ++n)
        { 
            if(ii[n] == IndNull) break;
            if(ii[n].m()==1) _index1.push_back(ii[n]); else { dim *= ii[n].m(); _indexn[++rn] = ii[n]; } 
        }
        allocate(dim);
        set_unique_Real();
    }

    explicit ITensor(const IndexVal& iv, Real fac = 1) : rn(0)
    { 
        if(iv.ind.m()==1) _index1.push_back(iv.ind); else { _indexn[++rn] = iv.ind; } 
        allocate(iv.ind.m());  
        p->v(iv.i) = fac; 
        set_unique_Real(); 
    }

    ITensor(const IndexVal& iv1, const IndexVal& iv2) : rn(0)
    { 
        if(iv1.ind.m()==1) _index1.push_back(iv1.ind); else { _indexn[++rn] = iv1.ind; } 
        if(iv2.ind.m()==1) _index1.push_back(iv2.ind); else { _indexn[++rn] = iv2.ind; } 
        allocate(iv1.ind.m()*iv2.ind.m());  
        p->v((iv2.i-1)*iv1.ind.m()+iv1.i) = 1; 
        set_unique_Real(); 
    }

    ITensor(IndexVal iv1, IndexVal iv2, IndexVal iv3,
            IndexVal iv4 = IVNull, IndexVal iv5 = IVNull, IndexVal iv6 = IVNull,
            IndexVal iv7 = IVNull, IndexVal iv8 = IVNull)
            : rn(0)
	{
        array<IndexVal,NMAX+1> iv = {{ IVNull, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        int dim = 1;
        for(int n = 1; n <= NMAX; ++n)
        { 
            if(iv[n].ind == IndNull) break;
            if(iv[n].ind.m()==1) _index1.push_back(iv[n].ind); else { dim *= iv[n].ind.m(); _indexn[++rn] = iv[n].ind; } 
        }
        allocate(dim);
        p->v((((((((iv8.i-1)*m(7)+iv7.i-1)*m(6)+iv6.i-1)*m(5)+iv5.i-1)*m(4)+iv4.i-1)*m(3)+iv3.i-1)*m(2)+iv2.i-1)*m(1)+iv1.i) = 1; 
        set_unique_Real();
    }

    explicit ITensor(const vector<Index>& I) : rn(0)
    {
        int alloc_size = 1;
        for(size_t n = 0; n < I.size(); ++n)
        {
            const Index& i = I[n];
            if(i == IndNull) Error("ITensor: null Index in constructor.");
            if(i.m()==1) _index1.push_back(i); 
            else 
            { 
                if(rn == NMAX) Error("ITensor(const vector<Index>& I): too many indices with m > 1");
                GET(_indexn,++rn) = i; 
                alloc_size *= i.m(); 
            }
        }
        allocate(alloc_size);
        set_unique_Real();
    }

    ITensor(const vector<Index>& I, const Vector& V) : p(new ITDat(V)), rn(0)
    {
        int alloc_size = fillFromIndices(I);
        if(alloc_size != V.Length()) 
        { Error("ITensor(vector<Index>,Vector): incompatible Index and Vector sizes"); }
        set_unique_Real();
    }

    ITensor(const vector<Index>& I, const ITensor& other) : p(other.p), rn(0), scale_(other.scale_)
    {
        int alloc_size = fillFromIndices(I);
        if(alloc_size != other.vec_size()) 
        { Error("ITensor(vector<Index>,ITensor): incompatible Index and ITensor sizes"); }
        set_unique_Real();
    }

    ITensor(const vector<Index>& I, const ITensor& other, Permutation P) : p(0), rn(0), scale_(other.scale_)
    {
        int alloc_size = fillFromIndices(I);
        if(alloc_size != other.vec_size()) 
        { Error("ITensor(vector<Index>,ITensor,Permutation): incompatible Index and ITensor sizes"); }
        if(P.is_trivial()) { p = other.p; }
        else 
        { 
            allocate(); 
            other.reshapeDat(P,p->v); 
        }
        set_unique_Real();
    }

    ITensor(ITmaker itm) : rn(1)
	{
        GET(_indexn,1) = IndReIm; allocate(2);
        if(itm == makeComplex_1)  { p->v(1) = 1; }
        if(itm == makeComplex_i)  { p->v(2) = 1; }
        if(itm == makeConjTensor) { p->v(1) = 1; p->v(2) = -1; }
        set_unique_Real();
	}

    ITensor(istream& s) { read(s); }

    //ITensor: Read/Write -------------------------------------------------------

    void read(istream& s)
    { 
        bool null_;
        s.read((char*) &null_,sizeof(null_));
        if(null_) { *this = ITensor(); return; }
        s.read((char*) &rn,sizeof(rn));
        scale_.read(s);
        size_t i1size = 0;
        s.read((char*) &i1size,sizeof(i1size));
        p = new ITDat(s);
        for(int j = 1; j <= rn; ++j) _indexn[j] = Index(s);
        _index1.reserve(i1size); for(size_t i = 1; i <= i1size; ++i) _index1.push_back(Index(s));
        set_unique_Real();
    }

    void write(ostream& s) const 
    { 
        bool null_ = is_null();
        s.write((char*) &null_,sizeof(null_));
        if(null_) return;
        s.write((char*) &rn,sizeof(rn));
        scale_.write(s);
        const size_t i1size = _index1.size();
        s.write((char*) &i1size,sizeof(i1size));
        p->write(s);
        for(int j = 1; j <= rn; ++j) _indexn[j].write(s);
        foreach(const Index& I, _index1) I.write(s);
    }

    //Operators -------------------------------------------------------

    ITensor& operator*=(const ITensor& other);
    ITensor operator*(ITensor other) const { other *= *this; return other; }

    ITensor& operator*=(const IndexVal& iv) { ITensor oth(iv); return operator*=(oth); } 
    ITensor operator*(const IndexVal& iv) const { ITensor res(*this); res *= iv; return res; }
    friend inline ITensor operator*(const IndexVal& iv, ITensor t) { return (t *= iv); }

    ITensor& operator*=(Real fac) { scale_ *= fac; return *this; }
    ITensor operator*(Real fac) const { ITensor res(*this); res *= fac; return res; }
    friend inline ITensor operator*(Real fac, ITensor t) { return (t *= fac); }

    ITensor& operator/=(Real fac) { scale_ /= fac; return *this; }
    ITensor operator/(Real fac) const { ITensor res(*this); res /= fac; return res; }
    friend inline ITensor operator/(Real fac, ITensor t) { return (t /= fac); }

    //operator/=(ITensor) is actually non-contracting product
    ITensor& operator/=(const ITensor& other);
    ITensor operator/(const ITensor& other) const { ITensor res(*this); res /= other; return res; }

    ITensor& operator+=(const ITensor& o);
    ITensor operator+(const ITensor& o) const { ITensor res(*this); res += o; return res; }

    ITensor& operator-=(const ITensor& o)
    {
        if(this == &o) { scale_ = 0; return *this; }
        scale_ *= -1; operator+=(o); scale_ *= -1; return *this; 
    }
    ITensor operator-(const ITensor& o) const { ITensor res(*this); res -= o; return res; }

    //Index Methods ---------------------------------------------------

    Index findtype(IndexType t) const
	{
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j).type() == t) return GET(_indexn,j);
        foreach(const Index& I,_index1)
        if(I.type() == t) return I;
        Error("ITensor::findtype failed."); return Index();
	}

    bool findtype(IndexType t,Index& I) const
	{
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j).type() == t)
        {
            I = GET(_indexn,j);
            return true;
        }
        foreach(const Index& J,_index1)
        if(J.type() == t)
        {
            I = J;
            return true;
        }
        return false;
	}

    int findindex(const Index& I) const
    {
        if(I.m() == 1) return findindex1(I);
        else           return findindexn(I);
        return 0;
    }

    int findindexn(const Index& I) const
	{
        if(I.m() == 1) return 0;
        for(int j = 1; j <= rn; ++j)
        if(GET(_indexn,j) == I) return j;
        return 0;
	}

    int findindex1(const Index& I) const
	{
        if(I.m() != 1) return 0;
        for(size_t j = 0; j < _index1.size(); ++j)
        if(_index1[j] == I) return j+rn+1;
        return 0;
	}

    //Checks that if this has an Index I, then it
    //also has I' (I' can have any primelevel)
    bool has_symmetric_nindices() const
    {
        for(int i = 1; i <= rn; ++i)
        for(int j = 1; j <= rn; ++j)
        {
            if(j == i) continue;
            if(GET(_indexn,i).noprime_equals(GET(_indexn,j))) break;
            if(j == rn) return false;
        }
        return true;
    }

    bool has_common_index(const ITensor& other) const
    {
        for(int j = 1; j <= rn; ++j)
        for(int k = 1; k <= other.rn; ++k)
        if(_indexn[j] == other._indexn[k]) return true;

        for(size_t j = 0; j < _index1.size(); ++j)
        for(size_t k = 0; k < other._index1.size(); ++k)
        if(_index1[j] == other._index1[k]) return true;

        return false;
    }
    
    bool hasindex(const Index& I) const
	{
        if(I.m() == 1) return hasindex1(I);
        else           return hasindexn(I);
        return false;
	}

    bool hasindexn(const Index& I) const
	{
        for(int j = 1; j <= rn; ++j)
        if(_indexn[j] == I) return true;
        return false;
	}

    bool hasindex1(const Index& I) const
	{
        for(size_t j = 0; j < _index1.size(); ++j)
        if(_index1[j] == I) return true;
        return false;
	}

    bool notin(const Index& I) const { return !(hasindexn(I) || hasindex1(I)); }

    template <class Iterable>
    void addindex1(const Iterable& indices) 
    { 
        if(_index1.empty()) { _index1.insert(_index1.begin(),indices.begin(),indices.end()); }
        else
        {
            _index1.insert(_index1.begin(),indices.begin(),indices.end());
            sort(_index1.begin(),_index1.end());
            vector<Index>::iterator it = unique(_index1.begin(),_index1.end());
            _index1.resize(it-_index1.begin());
        }
        set_unique_Real();
    }

    void addindex1(const Index& I) 
    { 
        if(I.m() != 1) { cerr << "I = " << I << "\n"; Error("ITensor::addindex1: m != 1"); }
        if(hasindex1(I)) return;
        _index1.push_back(I); 
        set_unique_Real();
    }

    void removeindex1(const Index& I) 
    { 
        vector<Index>::iterator it = find(_index1.begin(),_index1.end(),I);
        if(it == _index1.end()) Error("Couldn't find m == 1 Index to remove.");
        _index1.erase(it);
        set_unique_Real();
    }

    //Removes the jth index as found by findindex
    void removeindex1(int j) 
    { 
        vector<Index>::iterator it = _index1.begin()+(j-(rn+1));
        _index1.erase(it);
        set_unique_Real();
    }


    //Primelevel Methods ------------------------------------

    void noprime(PrimeType p = primeBoth)
	{
        for(int j = 1; j <= rn; ++j) _indexn[j].noprime(p);
        foreach(Index& I, _index1) I.noprime(p);
        set_unique_Real();
	}

    void doprime(PrimeType pt, int inc = 1)
	{
        for(int j = 1; j <= rn; ++j) _indexn[j].doprime(pt,inc);
        foreach(Index& I, _index1) I.doprime(pt,inc);
        set_unique_Real();
	}

    void primeall() { doprime(primeBoth,1); }
    void primesite() { doprime(primeSite,1); }
    void primelink() { doprime(primeLink,1); }

    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        for(int j = 1; j <= rn; ++j) GET(_indexn,j).mapprime(plevold,plevnew,pt);
        foreach(Index& I, _index1) I.mapprime(plevold,plevnew,pt);
        set_unique_Real();
	}

    void mapprimeind(const Index& I, int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        for(int j = 1; j <= rn; ++j) 
        if(I == _indexn[j])
        {
            _indexn[j].mapprime(plevold,plevnew,pt);
            set_unique_Real();
            return;
        }
        foreach(Index& J, _index1) 
        if(I == J)
        {
            J.mapprime(plevold,plevnew,pt);
            set_unique_Real();
            return;
        }
        Print(*this);
        Print(I);
        Error("ITensor::mapprimeind: index not found.");
	}

    void primeind(const Index& I, int inc = 1) { mapindex(I,I.primed(inc)); }

    void primeind(const Index& I, const Index& J)
	{ 
        mapindex(I,I.primed());
        mapindex(J,J.primed());
	}

    void noprimeind(const Index& I) { mapindex(I,I.deprimed()); }

    friend inline ITensor primed(ITensor A)
    { A.doprime(primeBoth,1); return A; }

    friend inline ITensor primesite(ITensor A)
    { A.doprime(primeSite,1); return A; }

    friend inline ITensor primelink(const ITensor& A)
    { ITensor res(A); res.doprime(primeLink,1); return res; }

    friend inline ITensor primeind(ITensor A, const Index& I)
    { A.mapindex(I,I.primed()); return A; }

    friend inline ITensor primeind(ITensor A, const Index& I1, const Index& I2)
    { A.mapindex(I1,I1.primed()); A.mapindex(I2,I2.primed()); return A; }

    friend inline ITensor deprimed(ITensor A)
    { A.noprime(); return A; }

    //Element Access Methods ----------------------------------------


    Real val0() const 
	{ assert(p != 0); assert(rn == 0); return p->v(1)*scale_; }
    /*
    Real& ncval0()
	{ 
        assert(p != 0); 
        assert(rn == 0);
        solo(); 
        return p->v(1); 
    }
    */

    Real val1(int i1) const
	{ assert(p != 0); assert(rn == 1); return p->v(i1)*scale_; }
    /*
    Real& ncval1(int i1)
	{ 
        assert(p != 0); 
        assert(rn == 1);
        solo(); 
        return p->v(i1); 
    }
    */

    Real val2(int i1,int i2) const
	{ assert(p != 0); assert(rn == 2); return p->v((i2-1)*m(1)+i1)*scale_; }
    /*
    Real& ncval2(int i1,int i2)
	{ 
        assert(p != 0); 
        assert(rn == 2);
        solo(); 
        return p->v((i2-1)*m(1)+i1); 
    }
    */

    Real val3(int i1,int i2,int i3) const
	{ assert(p != 0); assert(rn == 3); return p->v(((i3-1)*m(2)+i2-1)*m(1)+i1)*scale_; }
    /*
    Real& ncval3(int i1,int i2,int i3)
    { 
        assert(p != 0); 
        assert(rn == 3);
        solo(); 
        return p->v(((i3-1)*m(2)+i2-1)*m(1)+i1); 
    }
    */

    Real val4(int i1,int i2,int i3,int i4) const
	{ assert(p != 0); assert(rn == 4); return p->v((((i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1)*scale_; }

    Real val5(int i1,int i2,int i3,int i4,int i5) const
	{ assert(p != 0); assert(rn == 5); return p->v(((((i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1)*scale_; }

    Real val6(int i1,int i2,int i3,int i4,int i5,int i6) const
	{ assert(p != 0); assert(rn == 6); return p->v((((((i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1)*scale_; }

    Real val7(int i1,int i2,int i3,int i4,int i5,int i6,int i7) const
	{ assert(p != 0); assert(rn == 7); return p->v(((((((i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1)*scale_; }

    Real val8(int i1,int i2,int i3,int i4,int i5,int i6,int i7, int i8) const
	{ 
        assert(p != 0); 
        return p->v((((((((i8-1)*m(7)+i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1)*scale_; 
    }
    /*
    Real& ncval8(int i1,int i2,int i3,int i4,int i5,int i6,int i7, int i8) 
	{ 
        assert(p != 0); 
        assert(!_neg);
        assert(_logfac==0);
        solo(); 
        return p->v((((((((i8-1)*m(7)+i7-1)*m(6)+i6-1)*m(5)+i5-1)*m(4)+i4-1)*m(3)+i3-1)*m(2)+i2-1)*m(1)+i1); 
    }
    */

    Real& operator()(const IndexVal& iv1, const IndexVal& iv2 = IVNull, const IndexVal& iv3 = IVNull,
                     const IndexVal& iv4 = IVNull, const IndexVal& iv5 = IVNull, const IndexVal& iv6 = IVNull,
                     const IndexVal& iv7 = IVNull, const IndexVal& iv8 = IVNull)
	{
        array<IndexVal,NMAX+1> iv = {{ IVNull, iv1, iv2, iv3, iv4, iv5, iv6, iv7, iv8 }};
        vector<int> ja(NMAX+1,1);
        int numgot = 0;
        for(int k = 1; k <= rn; ++k) //loop over indices of this ITensor
        {
            for(int j = 1; j <= NMAX; ++j)  // loop over the given indices
            if(_indexn[k] == iv[j].ind) 
            { ++numgot; ja[k] = iv[j].i; break; }
        }
        if(numgot != rn) 
        {
            cerr << format("numgot = %d, rn = %d\n")%numgot%rn;
            Error("ITensor::operator(): Not enough indices");
        }
	    assert(p != 0); 
        solo(); 
        setScale(1);
        return p->v((((((((ja[8]-1)*m(7)+ja[7]-1)*m(6)+ja[6]-1)*m(5)+ja[5]-1)*m(4)+ja[4]-1)*m(3)+ja[3]-1)*m(2)+ja[2]-1)*m(1)+ja[1]);
	}

    //Methods for Mapping to Other Objects --------------------------------------

    // Assume *this and other have same indices but different order.
    // Copy other into *this, without changing the order of indices in either
    // operator= would put the order of other into *this
    void assignFrom(const ITensor& other)
    {
        if(this == &other) return;
        if(fabs(other.ur - ur) > 1E-12)
        {
            Print(*this); Print(other);
            Error("assignFrom: unique Real not the same"); 
        }
        Permutation P; getperm(other._indexn,P);
        scale_ = other.scale_;
        if(p->count() != 1) { p = new ITDat(); }
#ifdef DO_ALT
        else { p->alt.clear(); }
#endif
        other.reshapeDat(P,p->v);
    }

    void toMatrix11NoScale(const Index& i1, const Index& i2, Matrix& res) const;
    void toMatrix11(const Index& i1, const Index& i2, Matrix& res) const; //puts in lfac
    void fromMatrix11(const Index& i1, const Index& i2, const Matrix& res);

    // group i1,i2; i3,i4
    void toMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,Matrix& res) const;
    void fromMatrix22(const Index& i1, const Index& i2, const Index& i3, const Index& i4,const Matrix& res);

    // group i1,i2; i3
    void toMatrix21(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const;
    void fromMatrix21(const Index& i1, const Index& i2, const Index& i3, const Matrix& res);

    // group i1; i2,i3
    void toMatrix12(const Index& i1, const Index& i2, const Index& i3, Matrix& res) const;
    void fromMatrix12(const Index& i1, const Index& i2, const Index& i3, const Matrix& res);

    int vec_size() const { return p->v.Length(); }
    void assignToVec(VectorRef v) const
	{
        if(p->v.Length() != v.Length()) Error("ITensor::assignToVec bad size");
        v = p->v;
        v *= scale_;
	}
    void assignFromVec(const VectorRef& v)
	{
        if(p->v.Length() != v.Length()) Error("ITensor::assignToVec bad size");
        scale_ = 1;
        assert(p != 0);
        if(p->count() != 1) 
        { 
            intrusive_ptr<ITDat> new_p(new ITDat(v)); 
            p.swap(new_p); 
        }
        else
        {
            p->v = v;
#ifdef DO_ALT
            p->alt.clear();
#endif
        }
	}

    void reshapeDat(const Permutation& p, Vector& rdat) const;
    void reshapeTo(const Permutation& P, ITensor& res) const;

    void reshape(const Permutation& P)
    {
        if(P.is_trivial()) return;
        Vector newdat;
        this->reshapeDat(P,newdat);
        newdat *= scale_;
        assignFromVec(newdat);
    }


    //Other Methods -------------------------------------------------

    void Randomize() { solo(); p->v.Randomize(); }

    void SplitReIm(ITensor& re, ITensor& im) const
    {
        re = *this;
        im = *this;
        if(!is_complex()) { im *= 0; return; }
        re *= IndReIm(1);
        im *= IndReIm(2);
    }
    inline void conj() { if(!is_complex()) return; operator/=(ConjTensor); }

    inline bool is_zero() const { return (norm() < 1E-20); } 

    Real sumels() const { return p->v.sumels() * scale_; }

    Real norm() const { return Norm(p->v) * scale_; }

    void doNormLog()
	{
        Real f = Norm(p->v);
        if(fabs(f-1) < 1E-12) return;
        solo();
        if(f != 0) { p->v *= 1.0/f; scale_ *= f; }
	}

    void setScale(LogNumber newscale) const
	{
        if(scale_ == newscale) return;
        solo();
        //cerr << format("newscale = %.3E (logNum = %f, sign = %d)\n")%Real(newscale)%newscale.logNum()%newscale.sign();
        //LogNumber vscale = (scale_/newscale);
        //cerr << format("vscale = %.3E (logNum = %f, sign = %d)\n")%Real(vscale)%vscale.logNum()%vscale.sign();
        p->v *= (scale_/newscale);
        scale_ = newscale;
	}

    void print(string name = "",Printdat pdat = HideData) const 
    { printdat = (pdat==ShowData); cerr << "\n" << name << " =\n" << *this << "\n"; printdat = false; }

    friend inline ostream& operator<<(ostream & s, const ITensor & t)
    {
        s << "logscale(incl in elems) = " << t.scale().logNum() << ", r = " << t.r() << ": ";
        int i = 0;
        for(i = 1; i <= t.r_n(); ++i) { s << t.indexn(i) << (i != t.r_n() ? ", " : "; "); }
        foreach(const Index& I, t.index1()) { s << I << (i != t.r() ? ", " : ""); ++i; }
        if(t.is_null()) s << " (dat is null)\n";
        else 
        {
            s << format(" (L=%d,N=%.2f)\n") % t.vec_size() % t.norm();
            if(printdat)
            {
                const Real scale = t.scale_;
                const Vector& v = t.p->v;
                Counter c; t.initCounter(c);
                for(; c != Counter::done; ++c)
                {
                    Real val = v(c.ind)*scale;
                    if(fabs(val) > 1E-10)
                    { s << c << " " << val << "\n"; }
                }
            }
            else s << "\n";
        }
        return s;
    }

}; //ITensor


inline Real Dot(const ITensor& x, const ITensor& y, bool doconj = true)
{
    if(x.is_complex())
	{
        ITensor res = (doconj ? conj(x) : x); res *= y;
        if(res.r() != 1) Error("Bad Dot 234234");
        return res.val0();
	}
    else if(y.is_complex())
	{
        ITensor res = x; res *= y;
        if(res.r() != 1) Error("Bad Dot 37298789");
        return res.val0();
	}

    ITensor res = x; res *= y;
    if(res.r() != 0) 
	{ x.print("x"); y.print("y"); Error("bad Dot"); }
    return res.val0();
}

inline void Dot(const ITensor& x, const ITensor& y, Real& re, Real& im, bool doconj = true)
{
    if(x.is_complex())
	{
        ITensor res = (doconj ? conj(x) : x); res *= y;
        if(res.r() != 1) error("Bad Dot 334234");
        re = res(IndReIm(1));
        im = res(IndReIm(2));
        return;
	}
    else if(y.is_complex())
	{
        ITensor res = x; res *= y;
        if(res.r() != 1) error("Bad Dot 47298789");
        re = res(IndReIm(1));
        im = res(IndReIm(2));
        return;
	}
    if(x.r() != y.r()) 
	{
        cerr << "x = " << x << "\n";
        cerr << "y = " << y << "\n";
        Error("bad Dot 122414");
	}
    ITensor res = x; res *= y;
    if(res.r() != 0) 
	{
        cerr << "x = " << x << "\n";
        cerr << "y = " << y << "\n";
        Error("bad Dot 20234");
	}
    re = res.val0();
    im = 0;
}

inline ITensor operator*(const IndexVal& iv1, const IndexVal& iv2) { ITensor t(iv1); return (t *= iv2); }
inline ITensor operator*(const IndexVal& iv1, Real fac) { return ITensor(iv1,fac); }
inline ITensor operator*(Real fac, const IndexVal& iv) { return ITensor(iv,fac); }

// Given Tensors which represent operators (e.g. A(site-1',site-1), B(site-1',site-1), 
// Multiply them, fixing primes C(site-1',site-1)
template<class Tensor>
inline Tensor multSiteOps(Tensor a, Tensor b) // a * b  (a above b in diagram, unprimed = right index of matrix)
{
    a.mapprime(1,2,primeSite);
    a.mapprime(0,1,primeSite);
    Tensor res = a * b;
    res.mapprime(2,1,primeSite);
    return res;
}



#ifdef THIS_IS_MAIN
ITensor Complex_1(makeComplex_1), Complex_i(makeComplex_i), ConjTensor(makeConjTensor);
const Index& ITensor::ReImIndex = IndReIm;
#endif

#endif
