#ifndef __ITENSOR_ITENSOR_H
#define __ITENSOR_ITENSOR_H
#include "types.h"
#include "real.h"
#include "allocator.h"
#include "index.h"
#include "permutation.h"
#include ".profiling/prodstats.h"
#include ".profiling/count_copies.h"

using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::string;
using std::stringstream;
using std::pair;
using std::make_pair;

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
    boost::array<int,NMAX+1> n;
    boost::array<int,NMAX+1> i;
    int ind;

    Counter() : rn_(0)
	{
        n.assign(1); n[0] = 0;
        reset(0);
	}

    Counter(const boost::array<Index,NMAX+1>& ii,int rn,int r) { init(ii,rn,r); }

    void init(const boost::array<Index,NMAX+1>& ii, int rn, int r)
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
        //set 'done' condition
        if(i[rn_] > n[rn_]) reset(0);
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

    inline bool notDone() const { return i[1] != 0; }

    friend inline ostream& operator<<(ostream& s, const Counter& c)
    {
        s << "("; 
        for(int i = 1; i < c.r_; ++i)
            {s << c.i[i] << " ";} 
        s << c.i[c.r_] << ")";
        return s;
    }
};

//#define DO_ALT
#ifdef DO_ALT
struct PDat
{
    Permutation I; 
    Vector v;
    PDat(const Permutation& P_, const Vector& v_) 
    : I(P_.inverse()), v(v_) { }
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
    std::vector<PDat> alt;
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
    
    void print() const { std::cout << "ITDat: v = " << v; }

    inline void* operator new(size_t size) throw(std::bad_alloc)
        { return allocator.alloc(); }

    inline void operator delete(void* p) throw()
        { return allocator.dealloc(p); }

    friend class ITensor;
    ENABLE_INTRUSIVE_PTR(ITDat)
private:
    static DatAllocator<ITDat> allocator;
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
    typedef boost::array<Index,NMAX+1>::const_iterator index_it;
    static const Index& ReImIndex;
private:
    //mutable: const methods may want to reshape data
    mutable intrusive_ptr<ITDat> p; 
    int r_,rn_;
    //Indices, maximum of 8 (index_[0] not used), mutable to allow reordering
    mutable boost::array<Index,NMAX+1> index_; 
    Real ur;
    //mutable since e.g. scaleTo is logically const
    mutable LogNumber scale_; 

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
            p = new ITDat(*p);
            IF_COUNT_COPIES(++copycount;)
        }
	}
    
    void set_unique_Real()
	{
        ur = 0;
        for(int j = 1; j <= r_; ++j)
            { ur += index_[j].unique_Real(); }
	}

    void _construct1(const Index& i1)
    {
        assert(r_ == 1);
        if(i1.m() != 1) rn_ = 1; 
        index_[1] = i1;
        allocate(i1.m());
        set_unique_Real();
    }

    void _construct2(const Index& i1, const Index& i2)
    {
        assert(r_ == 2);
        if(i1.m()==1) 
        {
            index_[1] = i2; index_[2] = i1; 
            rn_ = (i2.m() == 1 ? 0 : 1);
        }
        else 
        { 
            index_[1] = i1; index_[2] = i2; 
            rn_ = (i2.m() == 1 ? 1 : 2); 
        }
        allocate(i1.m()*i2.m()); 
        set_unique_Real();
    }

    template<class IndexContainer>
    int fillFromIndices(const IndexContainer& I, int size)
    {
        r_ = size;
        assert(r_ <= NMAX);
        assert(rn_ == 0);
        int r1_ = 0;
        boost::array<const Index*,NMAX+1> index1_;
        int alloc_size = 1;
        for(int n = 0; n < r_; ++n)
        {
            const Index& i = I[n];
            DO_IF_DEBUG(if(i == IndNull) \
            Error("ITensor: null Index in constructor.");)
            if(i.m()==1) { GET(index1_,++r1_) = &i; }
            else         { GET(index_, ++rn_) = i; alloc_size *= i.m(); }
        }
        for(int l = 1; l <= r1_; ++l) index_[rn_+l] = *(index1_[l]);
        set_unique_Real();
        return alloc_size;
    }

    //Prefer to map via a Combiner
    //Though 'mapindex' is useful for tested, internal calls
    void mapindex(const Index& i1, const Index& i2)
	{
        assert(i1.m() == i2.m());
        for(int j = 1; j <= r_; ++j) 
        if(GET(index_,j) == i1) 
        {
            GET(index_,j) = i2;
            set_unique_Real();
            return;
        }
        Print(i1);
        Error("ITensor::mapindex: couldn't find i1.");
	}

    void getperm(const boost::array<Index,NMAX+1>& oth_index_, Permutation& P) const
    {
        for(int j = 1; j <= r_; ++j)
        {
            bool got_one = false;
            for(int k = 1; k <= r_; ++k)
            if(oth_index_[j] == index_[k])
            { P.from_to(j,k); got_one = true; break; }
            if(!got_one)
            {
                cerr << "j = " << j << "\n";
                Print(*this); cerr << "oth_index_ = \n";
                foreach(const Index& I, oth_index_) { cerr << I << "\n"; }
                Error("ITensor::getperm: no matching index");
            }
        }
    }

    friend void toMatrixProd(const ITensor& L, const ITensor& R, 
                             int& nsamen, int& cdim,
                             boost::array<bool,NMAX+1>& contractedL, 
                             boost::array<bool,NMAX+1>& contractedR, 
                             MatrixRefNoLink& lref, MatrixRefNoLink& rref);

    int _ind(int i1, int i2, int i3, int i4, 
             int i5, int i6, int i7, int i8) const;

    int _ind2(const IndexVal& iv1, const IndexVal& iv2) const;

    int _ind8(const IndexVal& iv1, const IndexVal& iv2, 
              const IndexVal& iv3, const IndexVal& iv4 = IVNull, 
              const IndexVal& iv5 = IVNull,const IndexVal& iv6 = IVNull,
              const IndexVal& iv7 = IVNull,const IndexVal& iv8 = IVNull)
    const;

public:

    //Accessor Methods ----------------------------------------------

    //depends on indices only, unordered:
    inline Real unique_Real() const { return ur; } 
    inline const Index& index(int j) const 
        { assert(j <= r_); return GET(index_,j); }
    inline int r() const { return r_; }
    inline int m(int j) const { assert(j <= r_); return GET(index_,j).m(); }

    inline bool is_null() const { return (p == 0); }
    inline bool is_not_null() const { return (p != 0); }
    inline bool is_complex() const { return hasindexn(IndReIm); }
    inline bool is_not_complex() const { return !hasindexn(IndReIm); }

    inline LogNumber scale() const { return scale_; }

    //Can be used for iteration over Indices in a foreach loop
    //e.g. foreach(const Index& I, t.index() ) { ... }
    inline const pair<index_it,index_it> index() const  
        { return make_pair(index_.begin()+1,index_.begin()+r_+1); }

    void initCounter(Counter& C) const { C.init(index_,rn_,r_); }

    //Constructors --------------------------------------------------

    ITensor() : p(0), r_(0), rn_(0), ur(0)  { }

    ITensor(Real val) : r_(0), rn_(0)
	{ 
        allocate(1);
        p->v = val;
        set_unique_Real();
    }

    explicit ITensor(const Index& i1) : r_(1), rn_(0)
	{ _construct1(i1); }

    ITensor(const Index& i1, Real val) : r_(1), rn_(0)
	{ _construct1(i1); p->v = val; }

    ITensor(const Index& i1, const Vector& V) 
    : p(new ITDat(V)), r_(1), rn_(0)
	{ 
        if(i1.m() != V.Length()) 
            Error("Mismatch of Index and Vector sizes.");
        if(i1.m() != 1) rn_ = 1;
        index_[1] = i1;
        set_unique_Real();
    }

    ITensor(Index i1,Index i2) : r_(2), rn_(0)
	{ _construct2(i1,i2); }

    //Create an ITensor as a matrix with 'a' on the diagonal
    ITensor(Index i1,Index i2,Real a) : r_(2), rn_(0)
    {
        _construct2(i1,i2);
        if(rn_ == 2) //then index order is i1, i2
        {
            const int nn = min(i1.m(),i2.m());
            for(int i = 1; i <= nn; ++i) 
                p->v((i-1)*i1.m()+i) = a;
        }
        else { p->v(1) = a; }
    }

    ITensor(Index i1,Index i2,const MatrixRef& M) : r_(2), rn_(0)
    {
        _construct2(i1,i2);
        if(i1.m() != M.Nrows() || i2.m() != M.Ncols()) 
            { Error("Mismatch of Index sizes and matrix."); }
        MatrixRef dref; 
        p->v.TreatAsMatrix(dref,i2.m(),i1.m()); 
        dref = M.t();
    }

    ITensor(Index i1, Index i2, Index i3,
            Index i4 = IndNull, Index i5 = IndNull, Index i6 = IndNull,
            Index i7 = IndNull, Index i8 = IndNull)
            : rn_(0)
    {
        boost::array<Index,NMAX> ii = {{ i1, i2, i3, i4, i5, i6, i7, i8 }};
        int size = 3;
        while(ii[size] != IndNull) ++size;
        int alloc_size = fillFromIndices(ii,size);
        allocate(alloc_size);
    }

    explicit ITensor(const IndexVal& iv, Real fac = 1) : r_(1), rn_(0)
    { 
        _construct1(iv.ind);
        p->v(iv.i) = fac; 
    }

    ITensor(const IndexVal& iv1, const IndexVal& iv2) : r_(2), rn_(0)
    { 
        _construct2(iv1.ind,iv2.ind);
        p->v((iv2.i-1)*iv1.ind.m()+iv1.i) = 1; 
    }

    ITensor(const IndexVal& iv1, const IndexVal& iv2, 
            const IndexVal& iv3, const IndexVal& iv4 = IVNull, 
            const IndexVal& iv5 = IVNull, const IndexVal& iv6 = IVNull, 
            const IndexVal& iv7 = IVNull, const IndexVal& iv8 = IVNull)
            : rn_(0)
	{
        //Construct ITensor
        boost::array<Index,NMAX+1> ii = 
            {{ iv1.ind, iv2.ind, iv3.ind, iv4.ind, iv5.ind, 
               iv6.ind, iv7.ind, iv8.ind }};
        int size = 3; while(size < NMAX && ii[size+1] != IVNull.ind) ++size;
        int alloc_size = fillFromIndices(ii,size);
        allocate(alloc_size);

        //Assign specified element to 1
        boost::array<int,NMAX+1> iv = 
            {{ iv1.i, iv2.i, iv3.i, iv4.i, iv5.i, iv6.i, iv7.i, iv8.i }};
        boost::array<int,NMAX+1> ja; ja.assign(1);
        for(int k = 1; k <= rn_; ++k) //loop over indices of this ITensor
        {
            for(int j = 0; j < size; ++j)  // loop over the given indices
            { if(index_[k] == ii[j]) { ja[k] = iv[j]; break; } }
        }
        p->v(_ind(ja[1],ja[2],ja[3],ja[4],ja[5],ja[6],ja[7],ja[8])) = 1;
    }

    explicit ITensor(const std::vector<Index>& I) : rn_(0)
    {
        int alloc_size = fillFromIndices(I,I.size());
        allocate(alloc_size);
    }

    ITensor(const std::vector<Index>& I, const Vector& V) 
    : p(new ITDat(V)), rn_(0)
    {
        int alloc_size = fillFromIndices(I,I.size());
        if(alloc_size != V.Length()) 
            { Error("incompatible Index and Vector sizes"); }
    }

    ITensor(const std::vector<Index>& I, const ITensor& other) 
    : p(other.p), rn_(0), scale_(other.scale_)
    {
        int alloc_size = fillFromIndices(I,I.size());
        if(alloc_size != other.vec_size()) 
            { Error("incompatible Index and ITensor sizes"); }
    }

    ITensor(const std::vector<Index>& I, const ITensor& other, Permutation P) 
    : p(0), rn_(0), scale_(other.scale_)
    {
        int alloc_size = fillFromIndices(I,I.size());
        if(alloc_size != other.vec_size()) 
            { Error("incompatible Index and ITensor sizes"); }
        if(P.is_trivial()) { p = other.p; }
        else               { allocate(); other.reshapeDat(P,p->v); }
    }

    ITensor(ITmaker itm) : r_(1), rn_(1)
	{
        GET(index_,1) = IndReIm; allocate(2);
        if(itm == makeComplex_1)  { p->v(1) = 1; }
        if(itm == makeComplex_i)  { p->v(2) = 1; }
        if(itm == makeConjTensor) { p->v(1) = 1; p->v(2) = -1; }
        set_unique_Real();
	}

    ITensor(istream& s) { read(s); }

    //ITensor: Read/Write ---------------------------------------------------

    void read(istream& s)
    { 
        bool null_;
        s.read((char*) &null_,sizeof(null_));
        if(null_) { *this = ITensor(); return; }
        s.read((char*) &r_,sizeof(r_));
        s.read((char*) &rn_,sizeof(rn_));
        for(int j = 1; j <= r_; ++j) index_[j].read(s);
        scale_.read(s);
        p = new ITDat(s);
        set_unique_Real();
    }

    void write(ostream& s) const 
    { 
        bool null_ = is_null();
        s.write((char*) &null_,sizeof(null_));
        if(null_) return;
        s.write((char*) &r_,sizeof(r_));
        s.write((char*) &rn_,sizeof(rn_));
        for(int j = 1; j <= r_; ++j) index_[j].write(s);
        scale_.write(s);
        p->write(s);
    }

    //Operators -------------------------------------------------------

    ITensor& operator*=(const ITensor& other);
    ITensor operator*(ITensor other) const { other *= *this; return other; }

    ITensor& operator*=(const IndexVal& iv) 
        { ITensor oth(iv); return operator*=(oth); } 
    ITensor operator*(const IndexVal& iv) const 
        { ITensor res(*this); res *= iv; return res; }
    friend inline ITensor operator*(const IndexVal& iv, ITensor t) 
        { return (t *= iv); }

    ITensor& operator*=(Real fac) { scale_ *= fac; return *this; }
    ITensor operator*(Real fac) const 
        { ITensor res(*this); res *= fac; return res; }
    friend inline ITensor operator*(Real fac, ITensor t) 
        { return (t *= fac); }

    ITensor& operator/=(Real fac) { scale_ /= fac; return *this; }
    ITensor operator/(Real fac) const 
        { ITensor res(*this); res /= fac; return res; }
    friend inline ITensor operator/(Real fac, ITensor t) 
        { return (t /= fac); }

    //operator/=(ITensor) is actually non-contracting product
    ITensor& operator/=(const ITensor& other);
    ITensor operator/(const ITensor& other) const 
        { ITensor res(*this); res /= other; return res; }

    ITensor& operator+=(const ITensor& o);
    ITensor operator+(const ITensor& o) const 
        { ITensor res(*this); res += o; return res; }

    ITensor& operator-=(const ITensor& o)
    {
        if(this == &o) { scale_ = 0; return *this; }
        scale_ *= -1; operator+=(o); scale_ *= -1; return *this; 
    }
    ITensor operator-(const ITensor& o) const 
        { ITensor res(*this); res -= o; return res; }

    //Index Methods ---------------------------------------------------

    Index findtype(IndexType t) const
	{
        for(int j = 1; j <= rn_; ++j)
        if(index_[j].type() == t) return index_[j];
        Error("ITensor::findtype failed."); return Index();
	}

    bool findtype(IndexType t, Index& I) const
	{
        for(int j = 1; j <= r_; ++j)
        if(index_[j].type() == t)
        {
            I = index_[j];
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
        for(int j = 1; j <= rn_; ++j)
        if(index_[j] == I) return j;
        return 0;
	}

    int findindex1(const Index& I) const
	{
        for(int j = rn_+1; j <= r_; ++j)
        if(index_[j] == I) return j;
        return 0;
	}

    bool has_common_index(const ITensor& other) const
    {
        for(int j = 1; j <= r_; ++j)
        for(int k = 1; k <= other.r_; ++k)
        if(index_[j] == other.index_[k]) return true;

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
        for(int j = 1; j <= rn_; ++j)
        if(index_[j] == I) return true;
        return false;
	}

    bool hasindex1(const Index& I) const
	{
        for(int j = rn_+1; j <= r_; ++j)
        if(index_[j] == I) return true;
        return false;
	}

    bool notin(const Index& I) const { return !hasindex(I); }

    template <class Iterable>
    void addindex1(const Iterable& indices) 
    { 
        assert((r_+(int)indices.size()) <= NMAX);
        for(size_t j = 0; j < indices.size(); ++j)
        { 
            assert(indices[j].m() == 1);
            assert(!hasindex1(indices[j]));
            index_[++r_] = indices[j]; 
        }
        set_unique_Real();
    }

    void addindex1(const Index& I) 
    { 
        assert(I.m() == 1);
        assert(r_ < NMAX);
        assert(!hasindex1(I));
        index_[++r_] = I;
        set_unique_Real();
    }

    //Removes the jth index as found by findindex
    void removeindex1(int j) 
    { 
        assert(j <= r_);
        assert(j > rn_);
        for(int k = j; k < r_; ++k) index_[k] = index_[k+1];
        --r_;
        set_unique_Real();
    }

    inline void removeindex1(const Index& I) 
        { removeindex1(findindex1(I)); }



    //Primelevel Methods ------------------------------------

    void noprime(PrimeType p = primeBoth)
	{
        for(int j = 1; j <= r_; ++j) index_[j].noprime(p);
        set_unique_Real();
	}

    void doprime(PrimeType pt, int inc = 1)
	{
        for(int j = 1; j <= r_; ++j) index_[j].doprime(pt,inc);
        set_unique_Real();
	}

    void primeall() { doprime(primeBoth,1); }
    void primesite() { doprime(primeSite,1); }
    void primelink() { doprime(primeLink,1); }

    void mapprime(int plevold, int plevnew, PrimeType pt = primeBoth)
	{
        for(int j = 1; j <= r_; ++j) index_[j].mapprime(plevold,plevnew,pt);
        set_unique_Real();
	}

    void mapprimeind(const Index& I, int plevold, int plevnew, 
                     PrimeType pt = primeBoth)
	{
        for(int j = (I.m() == 1 ? rn_+1 : 1); j <= r_; ++j) 
        if(index_[j] == I)
        {
            index_[j].mapprime(plevold,plevnew,pt);
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

    friend inline ITensor primeind(ITensor A, const Index& I1, 
                                   const Index& I2)
    { A.mapindex(I1,I1.primed()); A.mapindex(I2,I2.primed()); return A; }

    friend inline ITensor deprimed(ITensor A)
    { A.noprime(); return A; }

    //Element Access Methods ----------------------------------------

    Real val0() const 
	{ assert(p != 0); assert(rn_ == 0); return p->v(1)*scale_.real(); }

    Real val1(int i1) const
	{ assert(p != 0); assert(rn_ <= 1); return p->v(i1)*scale_.real(); }

    Real& operator()()
	{ 
        if(rn_ != 0)
        {
            cerr << boost::format("# given = 0, rn_ = %d\n")%rn_;
            Error("Not enough indices (requires all having m!=1)");
        }
        assert(p != 0); 
        solo(); 
        scaleTo(1);
        return p->v(1);
    }

    const Real operator()() const
	{ 
        if(rn_ != 0)
        {
            cerr << boost::format("# given = 0, rn_ = %d\n")%rn_;
            Error("Not enough indices (requires all having m!=1)");
        }
        assert(p != 0); 
        return scale_.real()*p->v(1);
    }

    Real& operator()(const IndexVal& iv1)
	{
        assert(r_ >= 1);
        if(rn_ > 1) 
        {
            cerr << boost::format("# given = 1, rn_ = %d\n")%rn_;
            Error("Not enough indices (requires all having m!=1)");
        }
	    assert(p != 0); 
        solo(); 
        scaleTo(1);
        return p->v(iv1.i);
	}

    const Real operator()(const IndexVal& iv1) const
	{
        assert(r_ >= 1);
        if(rn_ > 1) 
        {
            cerr << boost::format("# given = 1, rn_ = %d\n")%rn_;
            Error("Not enough indices (requires all having m!=1)");
        }
	    assert(p != 0); 
        return scale_.real()*p->v(iv1.i);
	}

    inline Real& operator()(const IndexVal& iv1, const IndexVal& iv2) 
    {
	    assert(p != 0); 
        solo(); 
        scaleTo(1);
        return p->v(_ind2(iv1,iv2));
    }

    inline const Real operator()(const IndexVal& iv1, 
                                 const IndexVal& iv2) const
	{
	    assert(p != 0); 
        return scale_.real()*p->v(_ind2(iv1,iv2));
    }

    inline Real& operator()(const IndexVal& iv1, const IndexVal& iv2, 
                    const IndexVal& iv3, const IndexVal& iv4 = IVNull, 
                    const IndexVal& iv5 = IVNull,const IndexVal& iv6 = IVNull,
                    const IndexVal& iv7 = IVNull,const IndexVal& iv8 = IVNull)
    {
	    assert(p != 0); 
        solo(); 
        scaleTo(1);
        return p->v(_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8));
    }

    inline const Real operator()(const IndexVal& iv1, const IndexVal& iv2, 
                    const IndexVal& iv3, const IndexVal& iv4 = IVNull, 
                    const IndexVal& iv5 = IVNull,const IndexVal& iv6 = IVNull,
                    const IndexVal& iv7 = IVNull,const IndexVal& iv8 = IVNull) const
	{
	    assert(p != 0); 
        return scale_.real()*p->v(_ind8(iv1,iv2,iv3,iv4,iv5,iv6,iv7,iv8));
    }

    //Methods for Mapping to Other Objects ----------------------------------

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
        Permutation P; getperm(other.index_,P);
        scale_ = other.scale_;
        if(p->count() != 1) { p = new ITDat(); }
#ifdef DO_ALT
        else { p->alt.clear(); }
#endif
        other.reshapeDat(P,p->v);
    }

    void groupIndices(const boost::array<Index,NMAX+1>& indices, int nind, 
                      const Index& grouped, ITensor& res) const;

    void expandIndex(const Index& small, const Index& big, 
                     int start, ITensor& res) const;

    void fromMatrix11(const Index& i1, const Index& i2, const Matrix& res);
    void toMatrix11NoScale(const Index& i1, const Index& i2, 
                           Matrix& res) const;
    void toMatrix11(const Index& i1, const Index& i2, Matrix& res) const;

    /*
    // group i1,i2; i3,i4
    void toMatrix22(const Index& i1, const Index& i2, 
                    const Index& i3, const Index& i4,Matrix& res) const;
    void fromMatrix22(const Index& i1, const Index& i2, 
                      const Index& i3, const Index& i4,const Matrix& res);

    // group i1,i2; i3
    void toMatrix21(const Index& i1, const Index& i2, 
                    const Index& i3, Matrix& res) const;
    void fromMatrix21(const Index& i1, const Index& i2, 
                      const Index& i3, const Matrix& res);

    // group i1; i2,i3
    void toMatrix12(const Index& i1, const Index& i2, 
                    const Index& i3, Matrix& res) const;
    void fromMatrix12(const Index& i1, const Index& i2, 
                      const Index& i3, const Matrix& res);
    */

    int vec_size() const { return p->v.Length(); }
    void assignToVec(VectorRef v) const
	{
        if(p->v.Length() != v.Length()) 
            Error("ITensor::assignToVec bad size");
        v = p->v;
        v *= scale_.real();
	}
    void assignFromVec(const VectorRef& v)
	{
        if(p->v.Length() != v.Length()) 
            Error("ITensor::assignToVec bad size");
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
        solo();
        Vector newdat;
        this->reshapeDat(P,newdat);
        p->v = newdat;
    }

    //Other Methods -------------------------------------------------

    void Randomize() { solo(); p->v.Randomize(); }

    void SplitReIm(ITensor& re, ITensor& im) const
	{
	re = *this; im = *this;
	if(!is_complex()) { im *= 0; return; }
	//re *= IndReIm(1); im *= IndReIm(2);

	re.mapindex(IndReIm,IndReImP);
	im.mapindex(IndReIm,IndReImP);
	re *= IndReImP(1);
	im *= IndReImP(2);
	}

    inline void conj() { if(!is_complex()) return; operator/=(ConjTensor); }

    inline bool is_zero() const { return (norm() < 1E-20); } 

    Real sumels() const { return p->v.sumels() * scale_.real(); }

    Real norm() const { return Norm(p->v) * scale_.real(); }

    void scaleOutNorm() const
	{
        Real f = Norm(p->v);
        if(fabs(f-1) < 1E-12) return;
        solo();
        if(f != 0) { p->v *= 1.0/f; scale_ *= f; }
	}

    void scaleTo(LogNumber newscale) const
	{
        if(scale_ == newscale) return;
        solo();
        if(newscale.isRealZero()) { p->v = 0; }
        else 
	    {
            scale_ /= newscale;
            p->v *= scale_.real();
	    }
        scale_ = newscale;
	}

    void print(string name = "",Printdat pdat = HideData) const 
    { 
        printdat = (pdat==ShowData); 
        cerr << "\n" << name << " =\n" << *this << "\n"; 
        printdat = false; 
    }

    friend ostream& operator<<(ostream & s, const ITensor & t);

}; //ITensor


Real Dot(const ITensor& x, const ITensor& y, bool doconj = true);

void Dot(const ITensor& x, const ITensor& y, Real& re, Real& im, 
                bool doconj = true);

inline ITensor operator*(const IndexVal& iv1, const IndexVal& iv2) 
    { ITensor t(iv1); return (t *= iv2); }
inline ITensor operator*(const IndexVal& iv1, Real fac) 
    { return ITensor(iv1,fac); }
inline ITensor operator*(Real fac, const IndexVal& iv) 
    { return ITensor(iv,fac); }

// Given Tensors which represent operators 
//(e.g. A(site-1',site-1), B(site-1',site-1), 
// Multiply them, fixing primes C(site-1',site-1)
// a * b  (a above b in diagram, unprimed = right index of matrix)
template<class Tensor>
inline Tensor multSiteOps(Tensor a, Tensor b) 
{
    a.mapprime(1,2,primeSite);
    a.mapprime(0,1,primeSite);
    Tensor res = a * b;
    res.mapprime(2,1,primeSite);
    return res;
}



#ifdef THIS_IS_MAIN
ITensor Complex_1(makeComplex_1), 
        Complex_i(makeComplex_i), 
        ConjTensor(makeConjTensor);
const Index& ITensor::ReImIndex = IndReIm;
#endif

#endif
